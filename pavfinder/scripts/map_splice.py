#!/usr/bin/env python
import argparse
import os
import sys
import datetime
import pysam
import pavfinder as pv
from pavfinder.transcriptome.exon_mapper import ExonMapper
from pavfinder.transcriptome.transcript import Transcript
from pavfinder.transcriptome.novel_splice_finder import extract_features, filter_events, corroborate_genome
from pavfinder.transcriptome.adjacency import Adjacency
from pavfinder.transcriptome.read_support import find_support

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="query to genome bam")
    parser.add_argument("query_fasta", type=str, help="path to query sequences")
    parser.add_argument("gtf", type=str, help="gtf file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("outdir", type=str, help="path to output directory")
    parser.add_argument("--debug", dest="debug", action="store_true", help="debug mode")
    parser.add_argument("--r2c", type=str, help="reads to genome bam")
    parser.add_argument("--nproc", type=int, help="number of processes. Default:4", default=4)
    parser.add_argument("--min_support", type=int, help="minimum read support. Default:4", default=4)
    parser.add_argument("--suppl_annot", type=str, nargs="+", help="supplementary annotation file(s) for checking novel splice events")
    parser.add_argument("--genome_bam", type=str, help="genome bam")
    parser.add_argument("--max_diff_splice", type=int, help="maximum number of base differences in splice motif. Default:1", default=1)
    
    args = parser.parse_args()
    return args

def create_pysam_bam(path):
    if path is not None and os.path.exists(path):
        return pysam.AlignmentFile(path)
    return None

def create_pysam_fasta(path):
    if path is not None and os.path.exists(path):
        return pysam.FastaFile(path)
    return None

def create_pysam_tabix(path):
    if path is not None and os.path.exists(path):
        return pysam.Tabixfile(path, parser=pysam.asGTF())
    return None

def main():
    args = parse_args()
    bam = create_pysam_bam(args.bam)
    query_fasta = create_pysam_fasta(args.query_fasta)
    genome_fasta = create_pysam_fasta(args.genome_fasta)
    transcripts_dict = Transcript.extract_transcripts(args.gtf)
    annot_tabix = create_pysam_tabix(args.gtf)
    genome_bam = None
    if args.genome_bam:
        genome_bam = create_pysam_bam(args.genome_bam)
    
    em = ExonMapper(annot_tabix,
                    transcripts_dict,
                    genome_fasta,
                    debug = args.debug)

    annots = [args.gtf]
    if args.suppl_annot:
        annots.extend(args.suppl_annot)
    accessory_known_features = extract_features(annots)
    
    mappings, junc_adjs, events = em.map_aligns(bam,
                                                query_fasta,
                                                genome_fasta,
                                                accessory_known_features=accessory_known_features,
                                                max_diff=args.max_diff_splice)
    juncs_merged = Adjacency.merge(junc_adjs)

    events_merged = Adjacency.merge(events)

    if args.r2c:
        all_adjs = []
        if juncs_merged:
            all_adjs.extend(juncs_merged)
        if events_merged:
            all_adjs.extend(events_merged)
        if all_adjs:
            find_support(all_adjs, args.r2c, args.query_fasta, num_procs=args.nproc, debug=args.debug)
        if events_merged:
            filter_events(events_merged, args.min_support)

    if genome_bam:
        corroborate_genome(events_merged, genome_bam)
    
    cmd = ' '.join(sys.argv)
    time = datetime.datetime.now().strftime("%Y-%m-%d:%H:%M:%S")
    software = '%s %s' % (pv.__name__, pv.__version__)
    em.output_mappings(mappings, '%s/mappings.tsv' % args.outdir)
    em.output_juncs(juncs_merged, '%s/junctions.bed' % args.outdir)    
    em.output_events(events_merged, '%s/novel_splicing.bedpe' % args.outdir, header=(software, '%s %s' % (time, cmd)))
    
main()
