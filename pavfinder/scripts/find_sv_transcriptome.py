#!/usr/bin/env python
import argparse
import os
import datetime
import pysam
import sys
from sets import Set
import pavfinder as pv
from pavfinder.transcriptome.transcript import Transcript
from pavfinder.transcriptome.sv_finder import SVFinder
from pavfinder.transcriptome.adjacency import Adjacency
from pavfinder.transcriptome.read_support import find_support

def combine_events(events, mappings):
    """Combine events via genome and transcripts alignment on contig level"""
    def same_event(event_t, event_g, window=100):
        if event_t.rearrangement == event_g.rearrangement or\
           event_t.event == event_g.event or\
           (event_t.event in ('fusion', 'read_through') and event_g.event in ('fusion', 'read_through')):
            seq_breaks_t = sorted(event_t.seq_breaks)
            if type(event_g.seq_breaks) is tuple or type(event_g.seq_breaks) is list:
                seq_breaks_g = sorted(event_g.seq_breaks)
            # after merging, tuple become str
            elif type(event_g.seq_breaks) is str:
                seq_breaks_g = sorted(map(int, event_g.seq_breaks.split('-')))

            if abs(seq_breaks_t[0] - seq_breaks_g[0]) <= window and\
               abs(seq_breaks_t[1] - seq_breaks_g[1]) <= window:
                if event_g.exon_bounds and\
                   event_g.exon_bounds[0] and\
                   event_g.exon_bounds[1]:
                    return event_g
                elif event_t.exon_bounds and\
                     event_t.exon_bounds[0] and\
                     event_t.exon_bounds[1]:
                    return event_t
                else:
                    if event_g.event in ('ITD') and not event_t.seq_id in contigs_same_genome_event:
                        return event_t
                    return event_g
        return False
            
    def same_mapping(query):
        passed = False
        if mappings['via_transcripts'] and mappings['via_genome']:
            if mappings['via_transcripts'].has_key(query) and\
               mappings['via_genome'].has_key(query):
                if not mappings['via_transcripts'][query] or not mappings['via_transcripts'][query][0]:
                    print '%s: mapping disagreed transcripts:None genome:%s' % (query,
                                                                                mappings['via_genome'][query])
                elif not mappings['via_genome'][query] or not mappings['via_genome'][query][0]:
                    print '%s: mapping disagreed transcripts:%s genome:None' % (query,
                                                                                mappings['via_transcripts'][query])
                elif mappings['via_transcripts'][query][0] and mappings['via_genome'][query][0]:
                    passed = True
                else:
                    print '%s: mapping disagreed transcripts:%s genome:%s' % (query,
                                                                              mappings['via_transcripts'][query],
                                                                              mappings['via_genome'][query])
            elif not mappings['via_transcripts'].has_key(query):
                print '%s: mapping disagreed transcripts:None genome:%s' % (query,
                                                                            mappings['via_genome'][query])
            else:
                print '%s: mapping disagreed transcripts:%s genome:None' % (query,
                                                                            mappings['via_transcripts'][query])
                
        return passed

    def extract_multiple_contig_events(events_by_contig):
        """Finds multiple contigs that map to same event, for use in same_event()"""
        events = []
        for contig in events_by_contig.keys():
            events.extend(events_by_contig[contig])
        events_merged = Adjacency.merge(events)
        contigs = Set()
        for event in events_merged:
            if ',' in event.seq_id:
                for contig in event.seq_id.split(','):
                    contigs.add(contig)
        return contigs

    combined_events = []
    contigs_same_genome_event = extract_multiple_contig_events(events['via_genome'])

    for query in Set(events['via_genome'].keys()) | Set(events['via_transcripts'].keys()):
        if events['via_genome'] and events['via_transcripts'] and not same_mapping(query):
            continue
        
        if not events['via_transcripts'].has_key(query):
            print 'onlygenome', query
            for event in events['via_genome'][query]:
                #if event.rearrangement == 'fusion' or event.rearrangement == 'read_through':
                combined_events.append(event)
            #else:
                #combined_events.extend(events['via_genome'][query])
        elif not events['via_genome'].has_key(query):
            print 'onlytranscripts', query
            for event in events['via_transcripts'][query]:
                combined_events.append(event)
        else:
            events_t = list(events['via_transcripts'][query])
            events_g = list(events['via_genome'][query])
            used_t = Set()
            used_g = Set()
            for i in range(len(events_t)):
                if i in used_t:
                    continue
                for j in range(len(events_g)):
                    if j in used_g:
                        continue
                    event = same_event(events_t[i], events_g[j])
                    if event:
                        combined_events.append(event)
                        used_t.add(i)
                        used_g.add(j)
                        
            for i in range(len(events_t)):
                if not i in used_t:
                    if not events_t[i].event in ('ins', 'del'):
                        combined_events.append(events_t[i])
            for i in range(len(events_g)):
                if not i in used_g:
                    if not events_g[i].event in ('ins', 'del'):
                        combined_events.append(events_g[i])
                        
    return combined_events
     
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("query_fasta", type=str, help="path to query sequences")
    parser.add_argument("gtf", type=str, help="gtf file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("outdir", type=str, help="path to output directory")
    parser.add_argument("--debug", dest="debug", action="store_true", help="debug mode")
    parser.add_argument("--tbam", type=str, help="query to transcript bam")
    parser.add_argument("--gbam", type=str, help="query to genome bam")
    parser.add_argument("--transcripts_fasta", type=str, help="path to target sequences")
    parser.add_argument("--r2c", type=str, help="reads to genome bam")
    parser.add_argument("--nproc", type=int, help="number of processes. Default:4", default=4)
    parser.add_argument("--genome_index", type=str, help="genome index path and name", nargs=2)
    parser.add_argument("--sort_by_coord", action="store_true", help="sort output by genome coordinates")
    parser.add_argument("--only_fusions", action="store_true", help="report only fusions and read-throughs")
    parser.add_argument("--version", action='version', version='%s %s' % (pv.__name__, pv.__version__))
    filtering = parser.add_argument_group('filtering')
    filtering.add_argument("--min_support", type=int, help="minimum read support. Default:4", default=4)
    filtering.add_argument("--min_overhang", type=int, help="minimum overhang for spanning reads. Default:4", default=4)
    filtering.add_argument("--min_indel_size", type=int, help="minimum indel size. Default:3", default=3)
    filtering.add_argument("--min_dup_size", type=int, help="minimum ins size to check for dup classification. Default:15", default=15)
    filtering.add_argument("--min_indel_flanking", type=int, help="minimum flanking contig lengths for indels. Default:10", default=10)
    filtering.add_argument("--no_utr", action="store_true", help="don't report events in UTR")
    filtering.add_argument("--include_nonsense_fusion", action="store_true", help="include non-sense fusions")
    filtering.add_argument("--include_non_exon_bound_fusion", action="store_true", help="include non-exon-bound fusions")
    filtering.add_argument("--include_noncoding_fusion", action="store_true", help="include noncoding fusions")
    filtering.add_argument("--max_homol_len", type=int, help="maximum homology sequence length. Default:10", default=10)
    filtering.add_argument("--max_novel_len", type=int, help="maximum novel sequence length. Default:10", default=10)
    filtering.add_argument("--subseq_len", type=int, help="subsequence length for filtering. Default:50", default=50)
    filtering.add_argument("--probe_len", type=int, help="probe sequence length for filtering. Default:100", default=100)
    filtering.add_argument("--disable_subseq_filtering", action="store_true", help="disable subseq filtering")
    filtering.add_argument("--disable_probe_filtering", action="store_true", help="disable probe filtering")

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
        
    gbam = create_pysam_bam(args.gbam)
    tbam = create_pysam_bam(args.tbam)
    query_fasta = create_pysam_fasta(args.query_fasta)
    genome_fasta = create_pysam_fasta(args.genome_fasta)
    transcripts_fasta = create_pysam_fasta(args.transcripts_fasta)
    transcripts_dict = Transcript.extract_transcripts(args.gtf)
    annot_tabix = create_pysam_tabix(args.gtf)
            
    sf = SVFinder(genome_fasta, annot_tabix, transcripts_dict, args.outdir, probe_len=args.probe_len, debug=args.debug)
    events = {'via_genome': {}, 'via_transcripts': {}}
    mappings = {'via_genome': {}, 'via_transcripts': {}}
    gene_hits = None
    if gbam and annot_tabix:
        events['via_genome'], mappings['via_genome'] = sf.find_events(gbam,
                                                                      query_fasta,
                                                                      genome_fasta,
                                                                      'genome',
                                                                      min_indel_size=args.min_indel_size,
                                                                      min_dup_size=args.min_dup_size,
                                                                      min_indel_flanking=args.min_indel_flanking,
                                                                      no_utr=args.no_utr,
                                                                      max_novel_len=args.max_novel_len,
                                                                      max_homol_len=args.max_homol_len,
                                                                      only_sense_fusion=not args.include_nonsense_fusion,
                                                                      only_exon_bound_fusion=not args.include_non_exon_bound_fusion,
                                                                      only_coding_fusion=not args.include_noncoding_fusion,
                                                                      only_fusions=args.only_fusions
                                                                      )
        
    if tbam:
        events['via_transcripts'], mappings['via_transcripts'] = sf.find_events(tbam,
                                                                                query_fasta,
                                                                                transcripts_fasta,
                                                                                'transcripts',
                                                                                external_mappings=mappings['via_genome'],
                                                                                min_indel_size=args.min_indel_size,
                                                                                min_dup_size=args.min_dup_size,
                                                                                min_indel_flanking=args.min_indel_flanking,
                                                                                no_utr=args.no_utr,
                                                                                no_indels=True,
                                                                                max_novel_len=args.max_novel_len,
                                                                                max_homol_len=args.max_homol_len,
                                                                                only_sense_fusion=not args.include_nonsense_fusion,
                                                                                only_exon_bound_fusion=not args.include_non_exon_bound_fusion,
                                                                                only_coding_fusion=not args.include_noncoding_fusion,
                                                                                only_fusions=args.only_fusions
                                                                                )        

    # combine events from genome and transcriptome alignments
    events_combined = combine_events(events, mappings)

    # merge identical events from different contigs
    events_merged = Adjacency.merge(events_combined)
    
    # filter by checking probe and subseq alignments
    if events_merged and args.genome_index and len(args.genome_index) == 2:
        if not args.disable_subseq_filtering:
            sf.filter_subseqs(events_merged, query_fasta, args.genome_index[0], args.genome_index[1], args.outdir,
                              subseq_len=args.subseq_len, debug=args.debug)
        if not args.disable_probe_filtering:
            sf.filter_probes(events_merged, args.genome_index[0], args.genome_index[1], args.outdir, args.probe_len, debug=args.debug)

    # read support
    if args.r2c:
        find_support(events_merged, args.r2c, args.query_fasta, min_overlap=args.min_overhang, num_procs=args.nproc, debug=args.debug)
        events_filtered = [event for event in events_merged if event.spanning >= args.min_support]
    else:
        events_filtered = events_merged

    # determine if events are in- or out-of-frame
    sf.set_frame(events_filtered, query_fasta, genome_fasta)

    # report (with meta data)
    cmd = ' '.join(sys.argv)
    time = datetime.datetime.now().strftime("%Y-%m-%d:%H:%M:%S")
    software = '%s %s' % (pv.__name__, pv.__version__)
    Adjacency.report_events(events_filtered, '%s/sv.bedpe' % args.outdir, sort_by_coord=args.sort_by_coord, header=(software, '%s %s' % (time, cmd)))

main()
    
