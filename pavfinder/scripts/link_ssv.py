#!/usr/bin/env python
from pavfinder.transcriptome.adjacency import Adjacency
import re
from pavfinder.transcriptome.novel_splice_finder import corroborate_genome, check_splice_motif_novel_site, check_splice_motif_skipped_exon
import pysam
import argparse
import sys
import datetime
from pavfinder.transcriptome.transcript import Transcript

def parse_splice(bedpe):
    headers = []
    events = []
    with open(bedpe, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                headers.append(line.rstrip())
            else:
                cols = line.rstrip().split('\t')
                if cols[0] == 'chrom1':
                    headers.append(line.rstrip())
                    fields = cols
                else:
                    event = {}
                    for field, col in zip(fields, cols):
                        event[field] = col
                    events.append(event)
    
    return headers, fields, events

def create_adjs(events, transcripts_dict, ref_fasta):
    p = re.compile('(chr\w+):(\d+)([A-Z]+)>([A-Z]+)')
    for event in events:
        adj = Adjacency(event['seq_id'], 
                        [event['chrom1'], event['chrom2']],
                        event['seq_breaks'], 
                        [event['end1'], event['start2']])
        adj.chroms = [event['chrom1'], event['chrom2']]
        adj.event = event['event']
        adj.exons = []
        for i in ('exon1', 'exon2'):
            if event[i].isnumeric():
                adj.exons.append(int(event[i]))
        adj.genome_support = None

        if event['splice_motif'] == 'na':
            adj.splice_motif = [None, []]
        else:
            adj.splice_motif = [event['splice_motif'], []]

        transcript = transcripts_dict[event['transcript']]

        if adj.event in ('novel_acceptor', 'novel_donor'):
            check_splice_motif_novel_site(adj, transcript, ref_fasta)
        elif adj.event == 'skipped_exon' and len(adj.exons) == 2:
            check_splice_motif_skipped_exon(adj, transcript, ref_fasta)

        if event['splice_site_variant'] != 'na':
            if ',' not in event['splice_site_variant']:
                m = p.match(event['splice_site_variant'])
            else:
                m = p.match(event['splice_site_variant'].split(',')[-1])

            if m:
                chrom = m.group(1)
                pos = m.group(2)
                from_base = m.group(3)
                to_base = m.group(4)

                adj.splice_motif = (event['splice_motif'], [(int(pos), to_base, event['splice_site_variant'].split(',')[-1])])

        event['adj'] = adj

def output(headers, fields, events, cmd, out_file):
    with open(out_file, 'w') as out:
        for header in headers:
            if header[0] == '#' and len(header.split()) > 2:
                out.write('{}\n'.format(cmd))
            else:
                out.write('{}\n'.format(header))
        for event in events:
            if 'adj' in event:
                if event['adj'].genome_support is not None:
                    if event['splice_site_variant'] == 'na':
                        event['splice_site_variant'] = event['adj'].splice_site_variant
                        event['genome_support'] = event['adj'].genome_support
                    else:
                        if not ',' in event['splice_site_variant']:
                            event['genome_support'] = event['adj'].genome_support
                        else:
                            event['genome_support'] = '{},{}'.format(event['adj'].genome_support, event['adj'].genome_support)
                    
            cols = []
            for field in fields:
                cols.append(event[field])
            out.write('{}\n'.format('\t'.join(cols)))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="input bedpe")
    parser.add_argument("output", type=str, help="output bedpe")
    parser.add_argument("gtf", type=str, help="gtf file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("--genome_bam", type=str, help="genome bam")
    parser.add_argument("--vcf", type=str, nargs="+", help="vcf(s)")
    parser.add_argument("--sample", type=str, help="sample in vcf file, assume first sample if not provided")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    headers, fields, events = parse_splice(args.input)

    transcripts_dict = Transcript.extract_transcripts(args.gtf)
    ref_fasta = pysam.FastaFile(args.genome_fasta)
    create_adjs(events, transcripts_dict, ref_fasta)
    adjs = [e['adj'] for e in events if 'adj' in e]

    genome_bam = None
    if args.genome_bam:
        genome_bam = pysam.AlignmentFile(args.genome_bam)
    
    vcfs = []
    if args.vcf:
        vcfs = [pysam.VariantFile(v) for v in args.vcf]

    corroborate_genome(adjs, genome_bam=genome_bam, vcfs=vcfs, sample=args.sample)

    cmd = ' '.join(sys.argv)
    time = datetime.datetime.now().strftime("%Y-%m-%d:%H:%M:%S")

    output(headers, fields, events, '#{} {}'.format(time, cmd), args.output)

main()
