import argparse
import re
import sys
import os
from collections import OrderedDict, defaultdict
from sets import Set
import subprocess
from pavfinder.shared.read_support import sum_support, filter_support

def parse_adjs(tsv):
    adjs = []
    with open(tsv, 'r') as ff:
        for line in ff:
            cols = line.rstrip('\n').split('\t')
            if cols[0] == 'ID':
                header = cols
            else:
                adj = OrderedDict()
                for i in range(len(cols)):
                    adj[header[i]] = cols[i]
                adjs.append(adj)
    return adjs

def parse_variants(vcf):
    variants = []
    with open(vcf, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                variants.append(line)
                if line[1:6] == 'CHROM':
                    header = line[1:].rstrip('\n').split('\t')
            else:
                variant = OrderedDict()
                cols = line.rstrip('\n').split('\t')
                for i in range(len(cols)):
                    if header[i] == 'INFO':
                        info = OrderedDict()
                        for pair in cols[i].split(';'):
                            attr, value = pair.split('=')
                            if attr in 'READSUPPORT':
                                continue
                            info[attr] = value
                        variant[header[i]] = info
                    else:
                        variant[header[i]] = cols[i]
                variants.append(variant)

    return variants

def create_coords_file(adjs, vid_to_aid, out_file):
    align_types = {}
    for vid, aids in vid_to_aid.iteritems():
        #print 'kk', vid, aids
        if len(aids) > 1:
            for aid in aids:
                align_types[aid] = 'split'
        else:
            align_types[aids[0]] = 'gapped'
    
    coords = defaultdict(list)
    for adj in adjs:
        for contig, start, end in zip(adj['contig'].split(','),
                                      adj['contig-break1'].split(','),
                                      adj['contig-break2'].split(',')):
            coords[contig].append((start, end, adj['ID']))
            
    with open(out_file, 'w') as out:
        for contig in sorted(coords.keys()):
            for start, end, aid in coords[contig]:
                #print 'gg', contig, start, end, aid
                out.write('%s\n' % '\t'.join([contig, start, end, align_types[aid]]))

def link_variants_adjs(variants, adjs):
    vid_to_aid = {}
    for v in range(len(variants)):
        if type(variants[v]) is str:
            continue
        aids = []
        for aid in variants[v]['ID'].split('-'):
            if aid[-1].isdigit():
                aids.append(aid)
            else:
                aids.append(aid[:-1])
            vid_to_aid[v] = aids
    return vid_to_aid

def parse_support(support_output):
    support = defaultdict(dict)
    with open(support_output, 'r') as ff:
        for line in ff:
            cols = line.rstrip().split('\t')
            spanning = int(cols[2])
            if cols[3].isdigit():
                flanking = int(cols[3])
            else:
                flanking = '-'
            support[cols[0]][cols[1]] = spanning, flanking
            
    return support

def filter_events(adjs, variants, vid_to_aid, support, min_support):
    failed_adjs = Set()
    adjs_by_id = {}
    for adj in adjs:
        adjs_by_id[adj['ID']] = adj
        spanning, flanking = [], []
        for contig, start, end in zip(adj['contig'].split(','),
                                      adj['contig-break1'].split(','),
                                      adj['contig-break2'].split(',')):
            span = '-'.join(map(str, sorted([start, end])))
            if support.has_key(contig) and support[contig].has_key(span):
                if type(support[contig][span][0]) is int:
                    spanning.append(support[contig][span][0])
                if type(support[contig][span][1]) is int:
                    flanking.append(support[contig][span][1])
        if spanning:
            adj['spanning_reads'] = sum_support(spanning)
        if flanking:
            adj['flanking_pairs'] = sum_support(flanking)
        else:
            adj['flanking_pairs'] = None

        if not filter_support(adj['spanning_reads'], adj['flanking_pairs'], min_support, use_spanning=True):
            failed_adjs.add(adj['ID'])

      
    failed_variants = Set()
    for i in range(len(variants)):
        if type(variants[i]) is str:
            continue

        aids = vid_to_aid[i]
        spanning, flanking = [], []
        
        skip = False
        for aid in aids:
            if not adjs_by_id.has_key(aid):
                skip = True
        if skip:
            continue
        
        for aid in aids:
            adj = adjs_by_id[aid]

            if aid in failed_adjs:
                failed_variants.add(i)
                break

            spanning.append(adj['spanning_reads'])
            flanking.append(adj['flanking_pairs'])

        if spanning:
            variants[i]['INFO']['SPANNING_READS'] = sum_support(spanning, use_minimum=True)
        if flanking:
            variants[i]['INFO']['FLANKING_PAIRS'] = sum_support(flanking, use_minimum=True)
        else:
            variants[i]['INFO']['FLANKING_PAIRS'] = None
        if not i in failed_variants:
            if not filter_support(variants[i]['INFO']['SPANNING_READS'], variants[i]['INFO']['FLANKING_PAIRS'], min_support, use_spanning=True):
                failed_variants.add(i)
                
    return failed_adjs, failed_variants

def get_support(coords_file, bam_file, contigs_fa, support_file, num_procs,
                allow_clipped_support=False,
                support_min_mapped=None):
    import pavfinder.shared.read_support
    script = '%s/read_support.py' % os.path.dirname(os.path.abspath(pavfinder.shared.read_support.__file__))
    cmd = "python %s %s %s %s %s -n %d" % (script,
                                           coords_file,
                                           bam_file,
                                           contigs_fa,
                                           support_file,
                                           num_procs)
    if allow_clipped_support:
        cmd += ' --allow_clipped_support'
    if support_min_mapped is not None:
        cmd += ' --support_min_mapped %s' % support_min_mapped
    print cmd
    #process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process = subprocess.Popen(cmd, shell=True)
    process.wait()
        
def output(adjs, variants, adjs_failed, variants_failed, outdir, out_prefix=None):
    if out_prefix is None:
        out_file = '%s/adjacencies_filtered.tsv' % outdir
    else:
        out_file = '%s/adjacencies_%s_filtered.tsv' % (outdir, out_prefix)
    
    with open(out_file, 'w') as tsv:
        header = adjs[0].keys()
        tsv.write('\t'.join(header) + '\n')

        by_id = {}
        for adj in adjs:
            if adj['ID'] in adjs_failed:
                continue

            cols = []
            for h in header:
                cols.append(str(adj[h]))
            tsv.write('\t'.join(cols) + '\n')
            by_id[cols[0]] = adj
        
    if out_prefix is None:
        out_file = '%s/variants_filtered.vcf' % outdir
    else:
        out_file = '%s/variants_%s_filtered.vcf' % (outdir, out_prefix)
    vheader = None
    iheader = None
    with open(out_file, 'w') as vcf:
        for i in range(len(variants)):
            if type(variants[i]) is str:
                vcf.write(variants[i])
            elif not i in variants_failed:
                if vheader is None:
                    vheader = variants[i].keys()
                    
                cols = []
                for vh in variants[i].keys():
                    if vh != 'INFO':
                        cols.append(variants[i][vh])
                    else:
                        info = []
                        for ih in variants[i][vh].keys():
                            info.append('%s=%s' % (ih, variants[i][vh][ih]))
                        cols.append(';'.join(info))
                line = '\t'.join(cols)
                vcf.write('%s\n' % line)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sv_dir", type=str, help="directory of sv output")
    parser.add_argument("r2c", type=str, help="r2c bam file")
    parser.add_argument("contigs_fa", type=str, help="contig fasta file")
    parser.add_argument("--out_prefix", type=str, help="prefix of output file", default=None)
    parser.add_argument("--num_procs", help="Number of processes. Default: 5", default=5, type=int)
    parser.add_argument("--min_support", help="minimum read support. Default:2", type=int, default=2)
    parser.add_argument("--allow_clipped_support", help="allow using clipped reads in gathering read support", action="store_true", default=False)
    parser.add_argument("--support_min_mapped", help="when clipped reads are allowed as read support, minimum ratio of read length mapped", type=float, default=None)
    
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    adjs = parse_adjs(args.sv_dir + '/adjacencies.tsv')
    variants = parse_variants(args.sv_dir + '/variants.vcf')
    print 'done parsing events'
    
    vid_to_aid = link_variants_adjs(variants, adjs)
    print 'done linking events'
    
    out_dir = args.sv_dir
    coords_file = '%s/coords.tsv' % out_dir
    create_coords_file(adjs, vid_to_aid, coords_file)
    print 'done creating coords file'
    
    # get support
    if args.out_prefix is not None:
        support_file = '%s/%s_support.tsv' % (out_dir, args.out_prefix)
    else:
        support_file = '%s/support.tsv' % out_dir
    if not os.path.exists(support_file):
        get_support(coords_file, args.r2c, args.contigs_fa, support_file, args.num_procs,
                    allow_clipped_support=args.allow_clipped_support,
                    support_min_mapped=args.support_min_mapped)
    
    if os.path.exists(support_file):
        print 'done getting support'
        # parse support
        support = parse_support(support_file)
        print 'done parsing support'
    
        # filter support
        adjs_failed, variants_failed = filter_events(adjs, variants, vid_to_aid, support, args.min_support)
        print 'done filtering support'

        # output
        output(adjs, variants, adjs_failed, variants_failed, args.sv_dir, out_prefix=args.out_prefix)
        print 'done making filtered output'
    
main()