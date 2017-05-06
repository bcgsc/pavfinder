#!/usr/bin/env python
import argparse
import pysam
import re
import sys
import subprocess
import os
from pavfinder.transcriptome.alignment import reverse_complement
from pavfinder.transcriptome.transcript import Transcript
from pavfinder.transcriptome.adjacency import Adjacency
from pavfinder.transcriptome.translate import is_inframe
from collections import defaultdict
from itertools import groupby

def find_unmapped_mates(bam, outdir):
    def write_seq(fa1, fa2, read):
        if read.is_read1:
            out = fa1
        else:
            out = fa2
            
        if read.is_reverse:
            seq = reverse_complement(read.query_sequence)
        else:
            seq = read.query_sequence
            
        out.write('>%s\n%s\n' % (read.query_name, seq))
    
    out_fas = ['%s/unmapped_%d.fa' % (outdir, i) for i in range(1,3)]
    fa1 = open(out_fas[0], 'w')
    fa2 = open(out_fas[1], 'w')

    prev_aln = None
    for aln in bam.fetch(until_eof=True):
	if not aln.is_unmapped and len(aln.cigartuples) > 1 and 5 in [t[0] for t in aln.cigartuples]:
	    continue

        if prev_aln is not None:
            if aln.is_unmapped and not aln.mate_is_unmapped and prev_aln.query_name == aln.query_name:
                write_seq(fa1, fa2, aln)
                write_seq(fa1, fa2, prev_aln)
                
            elif not aln.is_unmapped and aln.mate_is_unmapped and prev_aln.query_name == aln.query_name:
                write_seq(fa1, fa2, aln)
                write_seq(fa1, fa2, prev_aln)         
        prev_aln = aln
        
    fa1.close()
    fa2.close()
    return out_fas

def align_unmapped_mates(fastas, ref, nthreads, outdir):
    out_bam = '%s/unmapped.bam' % outdir
    cmd = 'bwa mem %s %s %s -t %s | samtools view -uhS - -o %s' % (ref,
                                                                   fastas[0],
                                                                   fastas[1],
                                                                   nthreads,
                                                                   out_bam)
    process = subprocess.Popen('/bin/bash -c "%s"' % cmd, shell = True)
    process.communicate()
    if process.returncode != 0:
	return Exception("Failed to run '%s'\n" % cmd)
    
    return out_bam

def find_breakpoints_from_clip(aln1, txt_seq2, orient2, txt_coord2):
    txt_break1, txt_break2 = None, None
    
    if aln1.cigartuples[0][0] == 4:
        clipped_seq = aln1.query_sequence[:aln1.cigartuples[0][1]]
        txt_break1 = aln1.reference_start + 1
    elif aln1.cigartuples[-1][0] == 4:
        clipped_seq = aln1.query_sequence[-1 * aln1.cigartuples[-1][1]:]
        txt_break1 = aln1.reference_end
    
    if orient2 == 'L':
        txt_seq_other_side = txt_seq2[txt_coord2:]
    elif orient2 == 'R':
        txt_seq_other_side = txt_seq2[:txt_coord2]
            
    matched_indices = [m.start() for m in re.finditer(clipped_seq, txt_seq_other_side)]
    if matched_indices:
        if orient2 == 'L':
            matched_index = matched_indices[0]
            distance = matched_index + 1
            txt_break2 = txt_coord2 + distance + len(clipped_seq)  - 1
    
        elif orient2 == 'R':
            matched_index = matched_indices[-1]
            distance = len(txt_seq_other_side) - matched_index
            txt_break2 = txt_coord2 - distance
    
    if txt_break1 is None or txt_break2 is None:
        return None
    
    return txt_break1, txt_break2

def set_genomic_orient(txt, txt_orient):
    if txt.strand == '+':
        orient = txt_orient
    else:
        if txt_orient == 'L':
            orient = 'R'
        else:
            orient = 'L'
    return orient
    
def shuffle_to_exon_bound(txt1, txt2, txt_break1, txt_break2, txt_orient1, txt_orient2, txt1_seq, txt2_seq, max_shift=10):
    def get_shift_distance(txt_break, orient, exon_span):
        if orient == 'L':
            return exon_span[1] - txt_break
        else:
            return txt_break - exon_span[0]
     
    exon_bound = False, False
    
    exon1 = txt1.txt_coord_to_exon(txt_break1)
    exon2 = txt2.txt_coord_to_exon(txt_break2)

    if exon1 is None or exon2 is None:
	return txt_break1, txt_break2, exon_bound

    exon_span1 = txt1.exon(exon1, transcript_coord=True)
    exon_span2 = txt2.exon(exon2, transcript_coord=True)
    shift1 = get_shift_distance(txt_break1, txt_orient1, exon_span1)
    shift2 = get_shift_distance(txt_break2, txt_orient2, exon_span2)
    
    if shift1 == 0 and shift2 == 0:
        exon_bound = True, True
    elif shift1 == 0:
        exon_bound = True, False
    elif shift2 == 0:
        exon_bound = False, True
    
    if exon_bound == (False, False):        
        if shift1 <= max_shift:
            txt_break2_shifted = txt_break2 + shift1
            exon2_shifted = txt2.txt_coord_to_exon(txt_break2_shifted)
            if exon2_shifted is None:
                return txt_break1, txt_break2, exon_bound
            exon_span2 = txt2.exon(exon2_shifted, transcript_coord=True)
            shift2_shifted = get_shift_distance(txt_break2_shifted, txt_orient2, exon_span2)          
            txt1_seq_shifted = txt1_seq[txt_break1 : txt_break1 + shift1].upper()
            txt2_seq_shifted = txt2_seq[txt_break2 - 1 : txt_break2 + shift1 - 1].upper()
            
            if txt1_seq_shifted == txt2_seq_shifted and shift2_shifted == 0:
                txt_break1 = txt_break1 + shift1
                txt_break2 = txt_break2_shifted
                exon_bound = True, True
                homol = txt1_seq_shifted

    if not exon_bound:
        if shift2 <= max_shift:
            txt_break1_shifted = txt_break1 - shift2
            exon1_shifted = txt1.txt_coord_to_exon(txt_break1_shifted)
            exon_span1 = txt1.exon(exon1_shifted, transcript_coord=True)
            shift1_shifted = get_shift_distance(txt_break1_shifted, txt_orient1, exon_span1)
            txt1_seq_shifted = txt1_seq[txt_break1 : txt_break1 + shift2].upper()
            txt2_seq_shifted = txt2_seq[txt_break2 - 1 : txt_break2 + shift2 - 1].upper()
            if txt1_seq_shifted == txt2_seq_shifted and shift1_shifted == 0:
                txt_break2 = txt_break2 - shift2
                txt_break1 = txt_break1_shifted
                exon_bound = True, True
                homol = txt1_seq_shifted
            
    return txt_break1, txt_break2, exon_bound

def check_frame(adj, genome_fasta):
    adj.in_frame = False
    chimeric_seq = adj.transcripts[0].get_sequence(genome_fasta)[:adj.transcript_breaks[0]] +\
                 adj.transcripts[1].get_sequence(genome_fasta)[adj.transcript_breaks[1] - 1:]
    if adj.transcripts[0].is_coding() and adj.transcripts[1].is_coding():
        chimeric_seq = adj.transcripts[0].get_sequence(genome_fasta)[:adj.transcript_breaks[0]] +\
                     adj.transcripts[1].get_sequence(genome_fasta)[adj.transcript_breaks[1] - 1:]
        in_frame = is_inframe(adj.transcripts[0],
                              adj.transcripts[1],
                              (adj.transcript_breaks[0], adj.transcript_breaks[0] + 1),
                              chimeric_seq,
                              genome_fasta)
        if type(in_frame) is tuple:
	    adj.in_frame = True
	else:
	    adj.in_frame = in_frame

def create_adj(event):
    seq_id = ','.join(event['spanning'] + event['flanking'])
    adj = Adjacency(seq_id,
                    (event['txts'][0].chrom, event['txts'][1].chrom),
                    ('na', 'na'),
                    event['genome_breaks'],
                    )
    adj.event = 'fusion'
    adj.sense_fusion = True
    adj.transcripts = event['txts']
    adj.upstream_transcript = event['txts'][0]
    adj.downstream_transcript = event['txts'][1]
    adj.orients = event['orients']
    adj.genome_breaks = event['genome_breaks']
    adj.chroms = (event['txts'][0].chrom, event['txts'][1].chrom)
    adj.transcript_breaks = event['txt_breaks']
    
    adj.exons_oriented = adj.exons = event['exons']
    adj.exon_bounds_oriented = event['exon_bounds']
    adj.exon_bounds = event['exon_bounds']
    adj.spanning = len(event['spanning'])
    adj.flanking = len(event['flanking'])
    adj.support = adj.spanning + adj.flanking
    
    return adj
    
def find_discordant_pairs(bam, transcripts_dict, genome_fasta, min_pairs=2):
    fusions = defaultdict(list)
    for aln in bam.fetch(until_eof=True):
        if not aln.is_unmapped and not aln.mate_is_unmapped and aln.reference_id != aln.next_reference_id:            
            ref1 = bam.getrname(aln.reference_id)
            ref2 = bam.getrname(aln.next_reference_id)
            txt1 = transcripts_dict[ref1]
            txt2 = transcripts_dict[ref2]
            gene1 = txt1.gene
            gene2 = txt2.gene
            
            if gene1 < gene2:
                genes = gene1, gene2
            else:
                genes = gene2, gene1
                
            fusions[genes].append(aln)
            
    adjs = []
    for genes in sorted(fusions, key=lambda p: len(fusions[p])):
        events = {}
        for read_id, pair in groupby(fusions[genes], lambda aln:aln.query_name):                
            pair_list = list(pair)
            alns1 = [aln for aln in pair_list if aln.is_read1 and not aln.is_unmapped]
            alns2 = [aln for aln in pair_list if not aln.is_read1 and not aln.is_unmapped]
            
            if len(alns1) == 1 and len(alns2) == 1:
                # 1==upstream 2==downstream
                aln1 = aln2 = None
                if alns1[0].is_reverse != alns2[0].is_reverse:
                    if not alns1[0].is_reverse:
                        aln1, aln2 = alns1[0], alns2[0]
                    else:
                        aln1, aln2 = alns2[0], alns1[0]

                if aln1 is None or aln2 is None:
                    continue
                
                ref1, ref2 = bam.getrname(aln1.reference_id), bam.getrname(aln2.reference_id)
                txt1, txt2 = transcripts_dict[ref1], transcripts_dict[ref2]
                txt_orient1, txt_orient2 = 'L', 'R'
                txt_break1, txt_break2 = aln1.reference_end, aln2.reference_start
                
                support = None
                exon_bound = False, False
                if aln1.query_length == aln1.query_alignment_length and\
                   aln2.query_length == aln2.query_alignment_length:
                    support = 'flanking'
                    
                else:
                    txt1_seq = txt1.get_sequence(genome_fasta)
                    txt2_seq = txt2.get_sequence(genome_fasta)
                    result = None
                    if aln1.cigartuples[0][0] == 4 or aln1.cigartuples[-1][0] == 4:
                        support = 'spanning'
                        result = find_breakpoints_from_clip(aln1, txt2_seq, txt_orient2, txt_break2)
                        if result is not None:
                            txt_break1, txt_break2 = result
                        
                    elif aln2.cigartuples[0][0] == 4 or aln2.cigartuples[-1][0] == 4:
                        support = 'spanning'
                        result = find_breakpoints_from_clip(aln2, txt1_seq, txt_orient1, txt_break1)
                        if result is not None:
                            txt_break2, txt_break1 = result
                        
                    if result is not None:
                        txt_break1_shifted, txt_break2_shifted, exon_bound = shuffle_to_exon_bound(txt1, txt2,
                                                                                                   txt_break1, txt_break2,
                                                                                                   txt_orient1, txt_orient2,
                                                                                                   txt1_seq, txt2_seq)
                        txt_break1 = txt_break1_shifted
                        txt_break2 = txt_break2_shifted
                        
                if support is not None:
                    orient1 = set_genomic_orient(txt1, txt_orient1)
                    orient2 = set_genomic_orient(txt1, txt_orient2)
                    exon1, exon2 = txt1.txt_coord_to_exon(txt_break1), txt2.txt_coord_to_exon(txt_break2)
                    if exon1 is None or exon2 is None:
                        continue
                    
                    gene1, gene2 = txt1.gene, txt2.gene
                    
                    if support == 'flanking':
                        txt_break1 = txt1.exon(exon1, transcript_coord=True)[1]
                        txt_break2 = txt2.exon(exon2, transcript_coord=True)[0]
                        
                    genome_break1, genome_break2 = txt1.txt_coord_to_genome_coord(txt_break1), txt2.txt_coord_to_genome_coord(txt_break2)
                    if genome_break1 is None or genome_break2 is None:
                        continue
                        
                    key = gene1, genome_break1, orient1, gene2, genome_break2, orient2
                    
                    if not events.has_key(key):
                        events[key] = {'txt_breaks': [txt_break1, txt_break2],
                                       'genome_breaks': [genome_break1, genome_break2],
                                       'orients': [orient1, orient2],
                                       'txts': [txt1, txt2],
                                       'exons': [exon1, exon2],
                                       'flanking':[],
                                       'spanning':[],
                                       'exon_bounds': (False, False),
                                       }
                    events[key][support].append(aln1.query_name)
                    if events[key]['exon_bounds'] != (True, True):
                        events[key]['exon_bounds'] = exon_bound                    

        num_pairs = 0
        for event in events.keys():
            num_pairs += len(events[event]['flanking'])
            num_pairs += len(events[event]['spanning'])
        
        if num_pairs >= min_pairs:
            for event in events.keys():
                adj = create_adj(events[event])
                check_frame(adj, genome_fasta)
                adjs.append(adj)
            
    return adjs

def cleanup(outdir):
    for ff in ('unmapped_1.fa', 'unmapped_2.fa', 'unmapped.bam'):
	os.remove('%s/%s' % (outdir, ff))

def parse_args():
    parser = argparse.ArgumentParser(description='Extracts discordant read pairs when fusion breakpoints are not captured')
    parser.add_argument("r2c", type=str, help="r2c_bam file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("gtf", type=str, help="gtf file")
    parser.add_argument("transcripts_fasta", type=str, help="path to target sequences")
    parser.add_argument("outdir", type=str, help="outdir")
    parser.add_argument("--min_pairs", type=int, help="minimum read pairs per fusion", default=2)
    parser.add_argument("--nthreads", type=int, help="number of threads for bwa mem", default=12)
    parser.add_argument("--no_cleanup", action='store_true', help="don't clean up intermediate fasta and bam file")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    # extract unmapped read pairs
    r2c_bam = pysam.AlignmentFile(args.r2c)
    unmapped_fastas = find_unmapped_mates(r2c_bam, args.outdir)
    
    # align unmapped reads
    unmapped_bam_file = align_unmapped_mates(unmapped_fastas, args.transcripts_fasta, args.nthreads, args.outdir)

    # extract fusions
    unmapped_bam = pysam.AlignmentFile(unmapped_bam_file)
    transcripts_dict = Transcript.extract_transcripts(args.gtf)
    adjs = find_discordant_pairs(unmapped_bam, transcripts_dict, pysam.FastaFile(args.genome_fasta), min_pairs=args.min_pairs)
    Adjacency.report_events(adjs, '%s/discordant_pairs.bedpe' % args.outdir)
    
    if not args.no_cleanup:
	cleanup(args.outdir)

main()

