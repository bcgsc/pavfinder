#!/usr/bin/env python
import pysam
import argparse
import subprocess
import os
from pavfinder.transcriptome.transcript import Transcript
from collections import defaultdict

def get_genes(list_file):
    genes = set()
    with open(list_file, 'r') as ff:
        for line in ff:
            genes.add(line.rstrip('\n'))
    return genes

def get_longest(transcripts, coding_first):
    sorted_by_size = sorted(transcripts, key = lambda transcript: transcript.length(), reverse=True)
    if coding_first:
        for transcript in sorted_by_size:
            if transcript.is_coding():
                return transcript

    return sorted_by_size[0]

def get_transcripts(annot, coding_only, genes=None, only_longest=False, coding_first=False):
    transcripts_dict = Transcript.extract_transcripts(annot)
    transcripts = []
    if only_longest:
        by_gene = defaultdict(list)
        for transcript in transcripts_dict.values():
            by_gene[transcript.gene].append(transcript)
        for gene in by_gene.keys():
            transcripts.append(get_longest(by_gene[gene], coding_first))
    else:
        transcripts = transcripts_dict.values()
    
    if genes and type(genes) is set:
        transcripts = [t for t in transcripts if t.gene in genes]
        captured_genes = set([t.gene for t in transcripts])
        
        for gene in genes - set([t.gene for t in transcripts]):
            print("can't find %s from gtf" % gene)
                
    if coding_only:
        transcripts = [t for t in transcripts if t.is_coding()]
        
    return transcripts

def get_genomic_sequence(by_gene, genome_fa):
    seq = {}
    for gene, transcripts in by_gene.items():
        start, end = None, None
        for transcript in transcripts:
            if start is None or transcript.exons[0][0] < start:
                start = transcript.exons[0][0]
            if end is None or transcript.exons[-1][1] > end:
                end = transcript.exons[-1][1]
                
        if start is not None and end is not None:
            seq[gene] = genome_fa.fetch(transcripts[0].chrom, start - 1, end)
            
    return seq
            
def output_single(out_file, transcripts, genome_fa, genomic_seq=None):
    out = open(out_file, 'w')
    for transcript in transcripts:
        out.write('>%s %s\n%s\n' % (transcript.id,
                                    transcript.gene,
                                    transcript.get_sequence(genome_fa).upper()
                                    )
                  )
        
    if genomic_seq is not None:
        for gene, seq in genomic_seq.items():
            out.write('>%s\n%s\n' % (gene, seq.upper()))

    out.close()
    
def output_by_gene(out_dir, by_gene, genome_fa, genomic_seq=None):
    for gene, transcripts in by_gene.items():
        with open('%s/%s.fa' % (out_dir, gene), 'w') as out:
            for transcript in transcripts:
                out.write('>%s %s\n%s\n' % (transcript.id,
                                            transcript.gene,
                                            transcript.get_sequence(genome_fa).upper()
                                            )
                          )
            if genomic_seq is not None and transcript.gene in genomic_seq:
                out.write('>%s\n%s\n' % (transcript.gene, genomic_seq[transcript.gene].upper()))
    
def bwa_index(fa):
    print('bwa', 'index', os.path.abspath(fa))
    subprocess.call(['bwa index %s' % fa], shell=True)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf", type=str, help="GTF file")
    parser.add_argument("out", type=str, help="output file for single output or directory name if output_by_gene")
    parser.add_argument("fa", type=str, help="indexed genome FASTA")
    parser.add_argument("--genes", type=str, help="text file of genes, one gene per line")
    parser.add_argument("--only_coding", action="store_true", help="only coding transcripts")
    parser.add_argument("--index", action="store_true", help="generate BWA index of output FASTA")
    parser.add_argument("--genomic", action="store_true", help="include genomic sequence")
    parser.add_argument("--output_by_gene", action="store_true", help="output by gene")
    parser.add_argument("--only_longest", action="store_true", help="only longest transcript")
    parser.add_argument("--coding_first", action="store_true", help="pick coding over non-coding transcript")
    
    args = parser.parse_args()
    
    return args

def main():
    args = parse_args()
    genes = None
    if args.genes:
        genes = get_genes(args.genes)
    transcripts = get_transcripts(args.gtf, args.only_coding, genes=genes, only_longest=args.only_longest, coding_first=args.coding_first)
    
    by_gene = defaultdict(list)
    if args.genomic or args.output_by_gene:
        for transcript in transcripts:
            by_gene[transcript.gene].append(transcript)
    
    fasta = pysam.FastaFile(args.fa)
    
    genomic_seq = None
    if args.genomic:
        genomic_seq = get_genomic_sequence(by_gene, fasta)
    
    if not args.output_by_gene:
        output_single(args.out, transcripts, fasta, genomic_seq=genomic_seq)
    
        if args.index:
            bwa_index(args.out)
            
    else:
        if os.path.isdir(args.out):
            out_dir = args.out
        else:
            out_dir = os.path.dirname(args.out)
        output_by_gene(out_dir, by_gene, fasta, genomic_seq=genomic_seq)

main()
    
    
