import re
import sys
import os
import glob
from sets import Set
from operator import itemgetter
import collections
from alignment import Alignment, reverse_complement, search_by_regex, search_by_align, has_canonical_target
from chimera import find_chimera
from adjacency import Adjacency
from translate import check_frame, nuc_to_aa
from exon_mapper import ExonMapper
from itertools import groupby
import subprocess
import pysam
from intspan import intspan

class SVFinder:
    
    def __init__(self, genome_fasta, annot, transcripts_dict, working_dir, probe_len=100, debug=False):
	self.genome_fasta = genome_fasta
	self.transcripts_dict = transcripts_dict
	self.annot = annot
	self.debug = debug
	self.exon_mapper = ExonMapper(annot, transcripts_dict, genome_fasta)
	self.working_dir = working_dir
	self.probe_len = probe_len
        
    def find_events(self, bam, query_fasta, target_fasta, target_type, 
                    gene_mappings=None, external_mappings=None,
                    min_indel_size=0, min_indel_flanking=0,
                    no_utr=False, no_indels=False, no_inv=True,
                    max_homol_len=5, max_novel_len=20,
                    only_sense_fusion=True, only_exon_bound_fusion=True,
                    only_coding_fusion=True,
                    ):
	def filter_adj(adj):
	    def check_junc_seq(homol=False, novel=False):
		junc_seq = None
		if homol:
		    junc_seq = adj.homol_seq
		    seq_type = 'homol'
		    max_len = max_homol_len
		elif novel:
		    junc_seq = adj.novel_seq
		    seq_type = 'novel'
		    max_len = max_novel_len
		    
		if junc_seq is not None:
		    if len(junc_seq) > max_len:
			print '%s: filter out %s - %s_seq %s too long' % (adj.seq_id,
			                                                  adj.event,
			                                                  seq_type,
			                                                  junc_seq)
			return False
	
		    if len(junc_seq) > 2 and self.is_homopolymer_fragment(junc_seq, min_pc=100):
			print '%s: filter out %s - %s_seq %s is homopolymer' % (adj.seq_id,
			                                                        adj.event,
			                                                        seq_type,
			                                                        junc_seq)
			return False
		return True
    
	    if no_indels and adj.rearrangement in ('ins', 'del'):
		return False
	    
	    if no_inv and adj.rearrangement in ('inv', 'inv-dup'):
		if adj.event is not None:
		    if adj.event != 'fusion':
			return False
		else:
		    return False

	    if no_utr and adj.feature:
		features = adj.feature.split(',')
		utrs = [feature for feature in features if 'UTR' in feature]
		if len(features) == len(utrs):
		    if self.debug:
			print '%s: %s filtered out - in UTR' % (adj.seq_id, adj.rearrangement)
		    return False
	    
	    if adj.transcripts and adj.exons and\
	       (adj.exons[0] is not None and adj.exons[1] is not None):
		if gene_mappings is not None:
		    if not screen_genes(adj):
			print '%s: event not mapped to same external gene' % adj.seq_id
			return False

	    if adj.transcripts is None or adj.transcripts[0] is None or adj.transcripts[1] is None:
		print '%s: cannot map event' % adj.seq_id
		return False
	    
	    if adj.exons is None or adj.exons[0] is None or adj.exons[1] is None:
		if not adj.event in ('fusion', 'read_through') or only_exon_bound_fusion:
		    print '%s: cannot map event to exon' % adj.seq_id
		    return False

	    if adj.homol_seq is not None and adj.homol_seq != 'na' and not check_junc_seq(homol=True):
		return False
	    
	    if adj.novel_seq is not None and adj.novel_seq != 'na' and not check_junc_seq(novel=True):
		return False
	    
	    if adj.probe is not None and len(adj.probe) > 0 and self.is_homopolymer_fragment(adj.probe):
		print '%s: probe sequence high in 1 base %s' % (adj.seq_id, adj.probe)
		return False
	    
	    if adj.event in ('fusion', 'read_through'):
		if adj.transcripts[0].gene == adj.transcripts[1].gene:
		    return False
		
		if adj.genome_breaks[0] >= adj.transcripts[0].txStart() and adj.genome_breaks[0] <= adj.transcripts[0].txEnd() and\
		   adj.genome_breaks[0] >= adj.transcripts[1].txStart() and adj.genome_breaks[0] <= adj.transcripts[1].txEnd():
		    return False
		
		if adj.genome_breaks[1] >= adj.transcripts[0].txStart() and adj.genome_breaks[1] <= adj.transcripts[0].txEnd() and\
		   adj.genome_breaks[1] >= adj.transcripts[1].txStart() and adj.genome_breaks[1] <= adj.transcripts[1].txEnd():
		    return False
		
		if only_coding_fusion and (not adj.transcripts[0].is_coding() or not adj.transcripts[1].is_coding()):
		    return False
		
		if (only_exon_bound_fusion or only_sense_fusion) and\
		   adj.sense_fusion is not None and not adj.sense_fusion:
		    print '%s: fusion not sense' % adj.seq_id
		    return False
		
		if only_exon_bound_fusion and adj.exon_bounds is not None and\
		   (type(adj.exon_bounds) is tuple or type(adj.exon_bounds) is list) and\
		   (not adj.exon_bounds[0] or not adj.exon_bounds[1]):
		    print '%s: not exon-bound fusions %s' % (adj.seq_id, adj.exon_bounds)
		    return False
		
	    return True

	def screen_genes(adj):
	    passed = [True, True]
	    for i in range(len(adj.transcripts)):
		if not adj.transcripts[i].gene in gene_mappings[adj.seq_id]:
		    passed[i] = False
	    if passed[0] == False and passed[0] == passed[1]:
		return False
	    else:
		return True
	    
	def process_split_aligns(aligns, query_seq, genes=None):
	    adjs = find_chimera(aligns, query_seq=query_seq, debug=self.debug)

	    adjs_filtered = []
	    for adj in adjs:
		# filtering out event if sequences overlap too much
		# otherwise is_fusion() will try to "adjust" exon_bounds to make it pass
		if type(adj.homol_seq) is str and len(adj.homol_seq) > max_homol_len:
		    print '%s: filter out %s - %s_seq %s too long' % (adj.seq_id,
		                                                      adj.rearrangement,
		                                                      adj.homol_seq,
		                                                      len(adj.homol_seq))
		    continue
		self.update_adj(adj, aligns, query_seq, target_type, block_matches=block_matches)
		if genes is not None:
		    add_genes_from_event(genes, adj)

		    if len(genes) > 1 and\
		       external_mappings is not None and\
		       external_mappings.has_key(adj.seq_id):
			if len(external_mappings[adj.seq_id][0]) == 1 and external_mappings[adj.seq_id][1] == 'full':
			    print '%s chimera mapped to single gene mapped:%s external:%s' % (adj.seq_id,
			                                                                      genes,
			                                                                      external_mappings[adj.seq_id][0]
			                                                                      )
			    continue

			if not genes.intersection(external_mappings[adj.seq_id][0]):
			    print '%s chimera not agreed mapped:%s external:%s' % (adj.seq_id,
			                                                           genes,
			                                                           external_mappings[adj.seq_id][0]
			                                                           )
			    continue

		if filter_adj(adj):
		    adjs_filtered.append(adj)

	    return adjs_filtered
	
	def compile_mapping(mapping, align):
	    if target_type == 'genome':
		chrom = align.target
		span = intspan('%s-%s' % (align.tstart, align.tend))
	    elif self.transcripts_dict.has_key(align.target):
		transcript = self.transcripts_dict[align.target]
		chrom = transcript.chrom
		span = intspan('%s-%s' % (transcript.exons[0][0], transcript.exons[-1][1]))
	    if not mapping.has_key(chrom):
		mapping[chrom] = []
	    mapping[chrom].append(span)
	    
	def add_genes_from_event(genes, event):
	    if event.transcripts:
		for transcript in event.transcripts:
		    if transcript is not None:
			genes.add(transcript.gene)
			
	def get_mapping_genes(block_matches, align):
	    tids = Set([tid for tid in block_matches.keys()])
	    full_matched_tids = [tid for tid in tids if self.exon_mapper.is_full_match(block_matches[tid], perfect=False)]
	    genes = Set([self.transcripts_dict[tid].gene for tid in tids])
	    genes_full = Set([self.transcripts_dict[tid].gene for tid in full_matched_tids])
	    if len(genes_full) > 1:
		tid = self.exon_mapper.pick_best_mapping(block_matches, align)
		genes = Set()
		genes.add(self.transcripts_dict[tid].gene)
		return [tid], genes
	    elif len(genes_full) == 1:
		return full_matched_tids, genes_full
	    else:
		return tids, genes

	events_by_query = {}
	mappings_by_query = {}
	partial_aligns = []
	for query, group in groupby(bam.fetch(until_eof=True), lambda aln: aln.query_name):
	    print 'processing', query
	    events = []
	    query_seq = query_fasta.fetch(query)
	    genes = Set()
	    tids = []

	    aligns = []
	    for aln in list(group):
		if not aln.is_unmapped:
		    align = Alignment.from_alignedRead(aln, bam)
		    # skip align if it's invalid (e.g. begin with softclip) or 
		    # query sequence is potentially homopolymer
		    if align.is_valid():
		    #if align.is_valid() and\
		       #align.aligned_seq(query_seq) and\
		       #not self.is_homopolymer_fragment(align.aligned_seq(query_seq)):
			if target_type == 'genome':
			    if align.has_canonical_target():
				aligns.append(align)
			else:
			    aligns.append(align)
		    
	    if not aligns:
		continue
	    
	    # partial
	    partially_aligned = None
	    if len(aligns) == 1:
		partially_aligned = aligns[0].is_partial(query_seq)
	    
	    aligns_mapped = 0
	    block_matches = None
	    if target_type == 'genome':
		if len(aligns) == 1:
		    block_matches = self.exon_mapper.map_align(aligns[0])
		
		    # cannot map alignment to gene, skip
		    if not block_matches:
			print '%s:cannot map to transcript' % aligns[0].query
			continue

		    tids, genes = get_mapping_genes(block_matches, aligns[0])
		    if tids and genes:
			aligns_mapped = 1
		    if len(genes) > 1:
			adj = self.is_read_through_from_single_align(block_matches, aligns[0])
			if adj is not None:
			    self.update_adj(adj, (aligns[0],aligns[0]), query_seq, target_type, block_matches=block_matches)
			    if filter_adj(adj):
				events.append(adj)
		    elif len(tids) == 1 and not partially_aligned and self.exon_mapper.is_full_match(block_matches[list(tids)[0]]):
			mappings_by_query[query] = genes, 'full'
			continue

		else:
		    for align in aligns:
			block_matches = self.exon_mapper.map_align(align)
			if block_matches:
			    tids_align, genes_align = get_mapping_genes(block_matches, align)
			    if len(genes_align) >= 1:
				aligns_mapped += 1
			    genes = genes.union(genes_align)
			# reset block matches for process_split_aligns()
			block_matches = None

	    if partially_aligned:
		if external_mappings is None or\
		   not external_mappings.has_key(query) or\
		   (external_mappings[query] and\
		    external_mappings[query][1] != 'full'):
		    partial_aligns.append((aligns[0], block_matches, partially_aligned))
    
	    # chimera
	    if len(aligns) > 1:
		aligns_sorted = sorted(aligns, key=lambda a: int(a.qstart))
		for align1, align2 in zip(aligns_sorted, aligns_sorted[1:]):
		    events.extend(process_split_aligns([align1, align2], query_seq, genes))
		        
	    # indels
	    for align in aligns:
		if target_type == 'transcripts':
		    if self.transcripts_dict.has_key(align.target):
			genes.add(self.transcripts_dict[align.target].gene)

		adjs = self.find_indels(align, query_fasta, target_fasta, target_type, 
		                        min_size=min_indel_size, min_flanking=min_indel_flanking, no_indels=no_indels)
		for adj in adjs:
		    self.update_adj(adj, (align, align), query_seq, target_type, block_matches=block_matches)
		    if target_type == 'genome' and not genes:
			add_genes_from_event(genes, adj)
		    
		    if filter_adj(adj):
			events.append(adj)

		    if adj.event in ('repeat_expansion', 'repeat_reduction'):
			if len(adj.repeat_seq) == 3 and\
			   adj.transcripts and\
			   not adj.transcripts[0].within_utr(adj.genome_breaks[0]) and\
			   not adj.transcripts[0].within_utr(adj.genome_breaks[1]):
			    adj.in_frame = True
			    self.adjust_for_amino_acid_repeat(adj)
			
	    if events:
		events_by_query[query] = events
		    
	    if genes:
		if not partially_aligned and aligns_mapped == len(aligns):
		    mappings_by_query[query] = genes, 'full'
		else:
		    mappings_by_query[query] = genes, 'partial'
		
	if partial_aligns:
	    new_aligns = self.map_partial_aligns(partial_aligns, query_fasta, target_type)
	    
	    for query, aligns in new_aligns.iteritems():
		query_seq = query_fasta.fetch(query)
		events = process_split_aligns(aligns, query_seq)
		if events:
		    events_by_query[query] = events

	return events_by_query, mappings_by_query
    
    @classmethod
    def filter_probes(cls, events, genome_index_dir, genome_index, working_dir, probe_length, debug=False):
	def create_query_fasta(events, fa_file, min_size=0):
	    qname_to_event = {}
	    fa = open(fa_file, 'w')
	    for event in events:
		# don't check splicing events
		if event.event is not None and\
		   event.event in ('alt_donor', 'alt_acceptor', 'skipped_exon'):
		    continue

		seq_ids = event.seq_id.split(',')
		probes = event.probe.split(',')
		for i in range(len(seq_ids)):
		    seq_id = seq_ids[i]
		    probe = probes[i]

		    if type(event.size) is str or event.size >= min_size:
			if len(probe) > 0:
			    qname = '%s:%s:%s' % (seq_id, event.key(), event.size)
			    fa.write('>%s\n%s\n' % (qname, probe))
			    qname_to_event[qname] = event
			else:
			    print 'probe empty', seq_id, event.keys(), probe
		
	    fa.close()
	    fai = fa_file + '.fai'
	    if os.path.exists(fai):
		os.remove(fai)
	    return qname_to_event
	
	def run_align(probes_fa, nthreads=12):
	    aln_bam_file = '%s/probes.bam' % working_dir
	    
	    cmd = 'gmap -D %s -d %s %s -n0 -f samse -t %d | samtools view -bhS - -o %s' % (genome_index_dir,
	                                                                                   genome_index,
	                                                                                   probes_fa,
	                                                                                   nthreads,
	                                                                                   aln_bam_file)
	    failed = False
	    try:
		if debug:
		    print cmd
		subprocess.call(cmd, shell=True)
	    except:
		sys.stderr.write('Failed to run:%s\n' % cmd)
		failed = True
		
	    failed = False
	    if not failed:
		return pysam.AlignmentFile(aln_bam_file)
	    else:
		return None
	    
	def parse_and_filter(bam, fasta_file, qname_to_event, indel_size_check=20):
	    def within_same_gene(alns, query_seq, event, min_mapped=0.9, window=50):
		"""check if subseq of fusion lie in same gene or do not lie in one of the genes"""
		query_spans = intspan()
		aligns = []
		for aln in alns:
		    if not aln.is_unmapped:
			align = Alignment.from_alignedRead(aln, bam)
			query_spans.add('%s-%s' % (align.qstart, align.qend))
			aligns.append(align)

		if float(len(query_spans)) / len(query_seq) > min_mapped:
		    mappings = collections.defaultdict(list)
		    for i in range(len(aligns)):
			for transcript in event.transcripts:
			    if aligns[i].target == transcript.chrom and\
			       ((aligns[i].tstart >= transcript.exons[0][0] - window and aligns[i].tstart <= transcript.exons[-1][1] + window)\
			        or\
			        (aligns[i].tend >= transcript.exons[0][0] - window and aligns[i].tend <= transcript.exons[-1][1] + window)):
				mappings[transcript.id].append(i)
		    if not mappings:
			return 'fusion subseq not mapped to either of the 2 genes'
		return False

	    def has_same_indel(aln, event_type, size):
		"""check if probe alignment contains same simple indel as expected"""
		if len(aln.cigartuples) >= 3:
		    for i in range(len(aln.cigartuples) - 2):
			if aln.cigartuples[i][0] == 0 and\
			   aln.cigartuples[i + 2][0] == 0 and\
			   ((event_type == 'ins' and aln.cigartuples[i + 1][0] == 1) or\
			    (event_type == 'del' and aln.cigartuples[i + 1][0] == 2)) and\
			   aln.cigartuples[i + 1][1] == size:
			    return True
		return False

	    fasta = pysam.FastaFile(fasta_file)
	    events_by_key = dict((events[i].key(), i) for i in range(len(events)))
    
	    remove = Set()
	    for query, group in groupby(bam.fetch(until_eof=True), lambda aln: aln.query_name):
		alns = list(group)
		event = qname_to_event[query]

		aln = alns[0]
		seq_id, key, size = aln.query_name.split(':')
		chrom = key.split('-')[1]
		event_type = key.split('-')[0]
		target = None
		if not aln.is_unmapped:
		    target = bam.getrname(aln.tid)
		bad = False		    
		failed_reason = ''
				
		if aln.is_unmapped:
		    failed_reason = 'cannot map probe'
		    bad = True

		elif event_type in ('ins', 'del'):
		    if aln.is_unmapped:
			failed_reason = 'indel probe not mapped'
			bad = True

		    elif size.isdigit() and int(size) <= indel_size_check:
			#if len(alns) > 1:
			    #failed_reason = 'chimeric probe aligns for simple indel'
			    #bad = True
			if not has_same_indel(aln, event_type, int(size)):
			    failed_reason = '%s %s probe align indel not matched %s' % (event_type, size, aln.cigarstring)
			    bad = True

		elif 'repeat' not in event_type: 
		    if event_type == 'fusion':
			if int(aln.get_tag('NM')) <= 2 and float(aln.query_alignment_length)/probe_length > 0.9:
			    failed_reason = 'probe aligned exclusively %d/%d=%.02f' % (aln.query_alignment_length,
			                                                               probe_length,
			                                                               float(aln.query_alignment_length)/probe_length)
			    bad = True

		    if not bad:
			if event_type in ('fusion', 'read_through'):
			    failed = within_same_gene(alns, fasta.fetch(aln.query_name), event)
			    if failed:
				failed_reason = failed
				bad = True
		    
		if bad:
		    print '%s probe failed: %s' % (query, failed_reason)
		    remove.add(events_by_key[key])
			
	    for i in sorted(list(remove), reverse=True):
		del events[i]
	    
	query_fa_file = '%s/probes.fa' % working_dir
	qname_to_event = create_query_fasta(events, query_fa_file)
	if qname_to_event:
	    bam = run_align(query_fa_file)
	    parse_and_filter(bam, query_fa_file, qname_to_event)
	    
	    if not debug:
		os.remove(bam.filename)
		for ff in glob.glob(query_fa_file + '*'):
		    os.remove(ff)

    @classmethod
    def filter_subseqs(cls, events, query_fa, genome_index_dir, genome_index, working_dir, subseq_len=None, debug=False):
	def create_query_fasta(events, fa_file, min_size=20):
	    count = 0
	    fa = open(fa_file, 'w')
	    for event in events:
		# don't check splicing events
		if event.event in ('alt_donor', 'alt_acceptor', 'skipped_exon'):
		    continue

		seq_id = event.seq_id.split(',')[0]
		seq_breaks = event.seq_breaks.split(',')[0].split('-')
		subseqs = event.get_subseqs(query_fa.fetch(seq_id), seq_breaks=seq_breaks, len_on_each_side=subseq_len)
		
		if 'repeat' in event.event:
		    continue
		
		for i in range(len(subseqs)):
		    fa.write('>%s:%s:%d\n%s\n' % (seq_id, event.key(), i, subseqs[i]))
		    count += 1
		
	    fa.close()
	    return count
	
	def run_align(probes_fa, nthreads=12):
	    aln_bam_file = '%s/subseqs.bam' % working_dir
	    
	    cmd = 'gmap -D %s -d %s %s -f samse -t %d | samtools view -bhS - -o %s' % (genome_index_dir,
	                                                                               genome_index,
	                                                                               probes_fa,
	                                                                               nthreads,
	                                                                               aln_bam_file)
	    failed = False
	    try:
		if debug:
		    print cmd
		subprocess.call(cmd, shell=True)
	    except:
		sys.stderr.write('Failed to run:%s\n' % cmd)
		failed = True
		
	    failed = False
	    if not failed:
		return pysam.AlignmentFile(aln_bam_file)
	    else:
		return None
	    
	def is_mapped(aln, min_aligned):
	    if not aln.is_unmapped and int(aln.get_tag('NM')) == 0:
		mapped_size = sum([t[1] for t in aln.cigartuples if t[0] == 0])
		if float(mapped_size)/aln.infer_query_length() >= min_aligned:
		    return True
		else:
		    return False

	def overlap(chrom1, span1, chrom2, pos2, window):
	    """ check 2 genomic positions overlap with a window """
	    if chrom1 == chrom2 and\
	       pos2 >= span1[0] - window and\
	       pos2 <= span1[1] + window:
		return True
	    return False
	
	def parse_and_filter(bam, window=100):
	    events_by_query = collections.defaultdict(list)
	    for i in range(len(events)):
		seq_id = events[i].seq_id.split(',')[0]
		events_by_query[seq_id].append(i)
    
	    remove = Set()
	    subseq_mappings = {}
	    subseq_align_tallies = {}
	    for query, group in groupby(bam.fetch(until_eof=True), lambda aln: aln.query_name):
		seq_id, key, part = query.split(':')
		event_type, chrom1, pos1, orient1, chrom2, pos2, orient2 = key.split('-')

		if event_type in ('ins', 'del'):
		    continue

		alns = list(group)
		
		mapped_alns = [aln for aln in alns if is_mapped(aln, 1.0)]
		if mapped_alns:
		    subseq_align_tallies[query] = 0
		    matched = []
		    if not subseq_mappings.has_key(seq_id):
			subseq_mappings[seq_id] = {}
		    if not subseq_mappings[seq_id].has_key(key):
			subseq_mappings[seq_id][key] = {}
		    for aln in mapped_alns:
			target = bam.getrname(aln.tid)

			if has_canonical_target(target):
			    subseq_align_tallies[query] += 1

			if overlap(target,
			           (aln.reference_start, aln.reference_end),
			           chrom1,
			           int(pos1),
			           window):
			    matched.append(1)

			if overlap(target,
			           (aln.reference_start, aln.reference_end),
			           chrom2,
			           int(pos2),
			           window):
			    matched.append(2)

		    subseq_mappings[seq_id][key][part] = matched

	    multimapped = [query for query in subseq_align_tallies if subseq_align_tallies[query] > 1]
			    
	    for seq_id in events_by_query.keys():
		if subseq_mappings.has_key(seq_id):
		    for i in events_by_query[seq_id]:
			# subseq will map to more than 1 region in dup events
			event_type = events[i].key().split('-')[0]
			if event_type in ('dup', 'ITD', 'PTD'):
			    continue

			# check for multimapping of both subseqs
			if '%s:%s:0' % (seq_id, events[i].key()) in multimapped or\
			   '%s:%s:1' % (seq_id, events[i].key()) in multimapped:
			    remove.add(i)
			    if '%s:%s:0' % (seq_id, events[i].key()) in multimapped:
				print '%s: remove %s - subseq 0 multimap' % (seq_id,
				                                             events[i].key())
			    elif '%s:%s:1' % (seq_id, events[i].key()) in multimapped:
				print '%s: remove %s - subseq 1 multimap' % (seq_id,
				                                             events[i].key())
			    break
			
	    for i in sorted(list(remove), reverse=True):
		del events[i]
	
	query_fa_file = '%s/subseqs.fa' % working_dir
	count = create_query_fasta(events, query_fa_file)
	if count > 0:
	    bam = run_align(query_fa_file)
	    parse_and_filter(bam)

	    if not debug:
		os.remove(bam.filename)
	if not debug:
	    os.remove(query_fa_file)
	
    def update_adj(self, adj, aligns, query_seq, target_type, block_matches=None):
	def fix_orients():
	    if target_type == 'transcripts':
		orients = []
		if adj.transcripts[0].strand == '-':
		    if adj.orients[0] == 'L':
			orients.append('R')
		    else:
			orients.append('L')
		else:
		    orients.append(adj.orients[0])

		if adj.transcripts[1].strand == '-':
		    if adj.orients[1] == 'L':
			orients.append('R')
		    else:
			orients.append('L')
		else:
		    orients.append(adj.orients[1])

		if orients:
		    adj.orients = tuple(orients)

	def sort_genome_breaks():
	    if not adj.is_genome_breaks_sorted():
		adj.reverse_genome_breaks(change=True, annotation=True)

	genes = None
	if target_type == 'genome':
	    if len(aligns) == 1:
		aligns = [aligns[0], aligns[0]]
	    adj.update_attrs(chroms = adj.targets,
                             genome_breaks = adj.target_breaks)
	    # dont want to overwrite for read_through
	    if not adj.transcripts:
		transcripts = self.map_transcripts_to_adj(adj, aligns, block_matches=block_matches)
		if transcripts and len(transcripts) == 2:
		    if transcripts[0] == transcripts[1]:
			genes = (transcripts[0].gene)
		    else:
			genes = (transcripts[0].gene, transcripts[1].gene)
	
		    adj.update_attrs(transcripts = transcripts)
		adj.update_transcript_breaks()
	    
	else:
	    # initial align order may be different from order in chimera
	    if len(aligns) > 1:
		adj_aligns = self.map_aligns_to_adj(adj, aligns)
	    else:
		adj_aligns = (aligns[0], aligns[0])
	    transcripts = [self.transcripts_dict[align.target] for align in adj_aligns]
	    if transcripts and len(transcripts) == 2:
		if transcripts[0] == transcripts[1]:
		    genes = (transcripts[0].gene)
		else:
		    genes = (transcripts[0].gene, transcripts[1].gene)
	
		adj.update_attrs(transcripts = transcripts,
		                 chroms = [transcript.chrom for transcript in transcripts],
		                 transcript_breaks = adj.target_breaks)
		adj.update_genome_breaks()
		
	if adj.transcripts:
	    adj.update_exons(target_type)
	    if genes and len(genes) > 1:
		self.is_fusion(adj, (aligns[0].strand, aligns[1].strand), target_type)
		    
	    self.update_feature(adj)
	    fix_orients()
	    sort_genome_breaks()
	    
	adj.set_probe(query_seq, len_on_each_side=self.probe_len/2)
	adj.size = adj.get_size()
	
	# if event is already defined, don't need to modify
	if adj.event is None:
	    if adj.rearrangement == 'del':
		self.is_del_alt_splicing_event(adj)

	    if adj.rearrangement == 'dup' and\
	       adj.transcripts and adj.exons is not None and\
	       adj.exons[0] is not None and adj.exons[1] is not None:
		if adj.transcripts[0].gene == adj.transcripts[1].gene and\
		   adj.transcripts[0].is_coding() and\
		   not adj.transcripts[0].within_utr(adj.genome_breaks[0]) and\
		   not adj.transcripts[0].within_utr(adj.genome_breaks[1]):
		    if adj.exon_bounds[0] and adj.exon_bounds[1]:
			adj.event = 'PTD'
		    else:
			adj.event = 'ITD'
			    
	if (adj.event == 'None' or adj.event is None) and adj.rearrangement is not None:
	    adj.event = adj.rearrangement

	# update support span
	if adj.event == 'ins' and adj.size > 20:
	    adj.support_span = (adj.seq_breaks[0], adj.seq_breaks[0] + 1)
	
	return genes
    
    def update_feature(self, adj):
	utrs = []
	for i in range(len(adj.genome_breaks)):
	    utr = adj.transcripts[i].within_utr(adj.genome_breaks[i])
	    if utr == 5:
		utrs.append('5UTR')
	    elif utr == 3:
		utrs.append('3UTR')
	if utrs:
	    if len(Set(utrs)) == 1:
		adj.feature = utrs[0]
	    elif utrs:
		adj.feature = ','.join(utrs)
	
    def map_aligns_to_adj(self, adj, aligns):
	def find_align_from_breakpoint(target, breakpoint, orient):
	    for i in range(len(aligns)):
		if target == aligns[i].target:
		    if orient == 'L' and breakpoint == aligns[i].tend:
			return i
		    if orient == 'R' and breakpoint == aligns[i].tstart:
			return i
	    return None
	
	aligns_mapped = []
	if aligns[0] == aligns[1]:
	    aligns_mapped = aligns
	else:
	    for i in range(2):
		align_idx = find_align_from_breakpoint(adj.targets[i],
		                                       adj.target_breaks[i],
		                                       adj.orients[i])
		if align_idx is not None:
		    aligns_mapped.append(aligns[align_idx])
	    
	return aligns_mapped
    
    def is_homopolymer_fragment(self, seq, min_pc=70):
	freq = collections.Counter(seq.upper())
	for base in ('A', 'G', 'T', 'C'):
	    if freq[base] * 100.0 / len(seq) >= min_pc:
		return True
	return False
    
    def map_transcripts_to_adj(self, adj, aligns, block_matches=None, genes=(None, None)):
	def compare_mapping(m1, m2):
	    if m1['exon_bound'] and not m2['exon_bound']:
		return -1
	    elif not m1['exon_bound'] and m2['exon_bound']:
		return 1
	    elif m1['coding'] and not m2['coding']:
		return -1
	    elif not m1['coding'] and m2['coding']:
		return 1
	    elif not m1['within_utr'] and m2['within_utr']:
		return -1
	    elif m1['within_utr'] and not m2['within_utr']:
		return 1
	    elif m1['num_matched_blocks'] > m2['num_matched_blocks']:
		return -1
	    elif m1['num_matched_blocks'] < m2['num_matched_blocks']:
		return 1
	    elif m1['num_perfect_matched_blocks'] > m2['num_perfect_matched_blocks']:
		return -1
	    elif m1['num_perfect_matched_blocks'] < m2['num_perfect_matched_blocks']:
		return 1
	    elif m1['num_exons'] > m2['num_exons']:
		return -1
	    elif m1['num_exons'] < m2['num_exons']:
		return 1
	    elif m1['len'] > m2['len']:
		return -1
	    elif m1['len'] < m2['len']:
		return 1
	    else:
		return 0

	def pick_best(chrom, breakpoint, orient, block_matches, gene):
	    metrics = {}
	    for tid, matches in block_matches.iteritems():
		transcript = self.transcripts_dict[tid]
		if gene is not None and transcript.gene != gene:
		    continue
		metrics[tid] = {'len': transcript.length(),
		                'num_exons': len(transcript.exons),
		                'coding': transcript.is_coding(),
		                'within_utr': transcript.within_utr(breakpoint),
		                }
		metrics[tid]['exon_bound'] = False
		if (orient == 'L' and matches[-1] is not None and matches[-1][-1][1][1] == '=') or\
		   (orient == 'R' and matches[0] is not None and matches[0][0][1][0] == '='):
		    metrics[tid]['exon_bound'] = True
		    
		metrics[tid]['num_matched_blocks'] = len([m for m in matches if m is not None and len(m) == 1])
		metrics[tid]['num_perfect_matched_blocks'] = len([m for m in matches if m is not None and
		                                                  len(m) == 1 and m[0][1] == '=='])
	
	    tids_sorted = sorted(metrics.keys(), cmp = lambda m1, m2: compare_mapping(metrics[m1], metrics[m2]))
		
	    return self.transcripts_dict[tids_sorted[0]]
	
	def map_transcript(align, chrom, genome_break, orient, block_matches=None, gene=None):
	    if block_matches is None:
		block_matches = self.exon_mapper.map_align(align)
	    transcript = None
	    if block_matches:
		transcript = pick_best(chrom,
		                       genome_break,
		                       orient,
		                       block_matches,
		                       gene)
	    return transcript
	
	transcripts = [None, None]
	aligns_mapped = self.map_aligns_to_adj(adj, aligns)
	if aligns_mapped and len(aligns_mapped) == len(aligns):
	    transcripts[0] = map_transcript(aligns_mapped[0],
	                                    adj.chroms[0],
	                                    adj.genome_breaks[0],
	                                    adj.orients[0],
	                                    block_matches = block_matches,
	                                    gene = genes[0])
	    if aligns_mapped[0] == aligns_mapped[1] and adj.event != 'read_through':
		transcripts[1] = transcripts[0]
	    else:
		transcripts[1] = map_transcript(aligns_mapped[1],
		                                adj.chroms[1],
		                                adj.genome_breaks[1],
		                                adj.orients[1],
		                                block_matches = block_matches,
		                                gene = genes[1],
		                                )			
	else:
	    print 'cannot map aligns to adj:%s' % adj.seq_id
		    
	if transcripts[0] is not None and transcripts[1] is not None:
	    return tuple(transcripts)
	return None
    
    def is_del_alt_splicing_event(self, event):
	if event.exon_bounds and event.transcripts:
	    if (event.exon_bounds[0] is False and event.exon_bounds[1] is True) or\
	       (event.exon_bounds[0] is True and event.exon_bounds[1] is False):
		if event.exon_bounds[0] is False:
		    genome_pos = event.genome_breaks[0]
		    unmatched_exon, matched_exon = event.exons[0], event.exons[1]
		else:
		    genome_pos = event.genome_breaks[1]
		    unmatched_exon, matched_exon = event.exons[1], event.exons[0]
		    
		if matched_exon < unmatched_exon:
		    motif = 'acceptor'
		else:
		    motif = 'donor'
			
		if motif == 'acceptor':
		    if event.transcripts[0].strand == '+':
			start, end = genome_pos - 2 - 1, genome_pos - 1
		    else:
			start, end = genome_pos, genome_pos + 2
		else:
		    if event.transcripts[0].strand == '+':
			start, end = genome_pos, genome_pos + 2
		    else:
			start, end = genome_pos - 2 - 1, genome_pos - 1
			
		bases = self.genome_fasta.fetch(event.transcripts[0].chrom, start, end)
		if event.transcripts[0].strand == '-':
		    bases = reverse_complement(bases)
		
		if motif == 'donor' and bases.upper() == 'GT':
		    event.event = 'alt_donor'
		elif motif == 'acceptor' and bases.upper() == 'AC':
		    event.event = 'alt_acceptor'
		    
	    elif event.exon_bounds[0] is True and event.exon_bounds[1] is True:
		event.event = 'skipped_exon'
    
    def is_homopolymer(self, seq):
	return len(Set(map(''.join, zip(*[iter(seq)] * 1)))) == 1
    
    def map_partial_aligns(self, partial_aligns, query_fasta, target_type, min_len=15):
	"""attempt to map partial aligns"""
	def extract_clipped_seq(query_seq, query_blocks):
	    """Extract sequences that are clipped at the end"""
	    clipped = []
	    # positive
	    if query_blocks[0][0] < query_blocks[0][1]:
		if query_blocks[0][0] > 1:
		    seq = query_seq[:query_blocks[0][0] - 1]
		    clipped.append((seq, 'start'))
		if query_blocks[-1][1] < len(query_seq):
		    seq = query_seq[query_blocks[-1][1]:]
		    clipped.append((seq, 'end'))
	    # negative
	    else:
		if query_blocks[-1][1] > 1:
		    seq = query_seq[:query_blocks[-1][1] - 1]
		    clipped.append((seq, 'start'))
		if query_blocks[0][0] < len(query_seq):
		    seq = query_seq[query_blocks[0][0]:]
		    clipped.append((seq, 'end'))
	    return clipped
	
	def create_alignment(query, query_span, target, target_span, strand):
	    align = Alignment(query = query,
		              target = target,
		              qstart = min(query_span),
		              qend = max(query_span),
		              tstart = target_span[0],
		              tend = target_span[1],
		              strand = strand)
	    align.query_len = align.qend - align.qstart + 1
	    align.blocks = [target_span]
	    align.query_blocks = [query_span]
	    return align

	def find_new_align(query, query_span, target, match, transcript):
	    if match[0] < match[1]:
		match_strand = '+'
	    else:
		match_strand = '-'
		    
	    target_span = sorted([match[0], match[1]])
	    if target_type == 'transcripts':
		strand = match_strand
	    else:
		if transcript.strand == '+':
		    strand = match_strand
		else:
		    if match_strand == '+':
			strand = '-'
		    else:
			strand = '+'
	    if target_type == 'genome':
		target_span = sorted([transcript.txt_coord_to_genome_coord(match[0]),
	                              transcript.txt_coord_to_genome_coord(match[1])])
		target = transcript.chrom
	    return create_alignment(query, query_span, target, target_span, strand)
	
	def get_transcript(block_matches):
	    transcript = None
	    if target_type == 'genome':
		if block_matches:
		    tid = self.exon_mapper.pick_best_mapping(block_matches, align)
		    if tid and self.transcripts_dict.has_key(tid):
			transcript = self.transcripts_dict[tid]
	    else:
		if self.transcripts_dict.has_key(align.target):
		    transcript = self.transcripts_dict[align.target]
	    return transcript
	
	def parse_partial_aligns(bam, query_target, query_span):
	    new_aligns = {}
	    for query_pos, group in groupby(bam.fetch(until_eof=True), lambda aln: aln.query_name):
		alns = list(group)
		matches = []
		query_len = None
		query, pos = query_pos.rsplit('-', 1)
		for aln in alns:
		    if aln.is_unmapped or aln.cigartuples[0][0] != 0 or aln.cigartuples[-1][0] != 0:
			continue
		    align = Alignment.from_alignedRead(aln, bam)
		    query_len = align.query_len
		    if align.is_valid and align.target == query_target[query]:
			tcoords_sorted = sorted([align.tstart, align.tend])
			if align.strand == '+':
			    matches.append(tcoords_sorted)
			else:
			    matches.append(tcoords_sorted[::-1])
		
		if matches and len(matches) == 1:
		    transcript = self.transcripts_dict[align.target]
		    match = matches[0]
		    if target_type == 'genome':
			target = transcript.chrom
		    else:
			target = align.target
		    align = find_new_align(query, query_span[query], target, match, transcript)
		    if align is not None:
			new_aligns[query] = align
	    return new_aligns
	
	def run_align(jobs, query_fa_file, target_fa_file):
	    query_fa = open(query_fa_file, 'w')
	    target_fa = open(target_fa_file, 'w')
	    target_seqs = {}
	    for align, qid, tid, query_seq, target_seq in jobs:
		target_seqs[tid] = target_seq
		query_fa.write('>%s\n%s\n' % (qid, query_seq))
		
	    for tid, seq in target_seqs.iteritems():
		target_fa.write('>%s\n%s\n' % (tid, seq))
		
	    query_fa.close()
	    target_fa.close()
	    
	    aln_bam_file = '%s/partial_aln.bam' % self.working_dir
	    index_cmd = 'bwa index %s' % target_fa_file
	    aln_cmd = 'bwa mem %s %s | samtools view -bhS - -o %s' % (target_fa_file,
	                                                              query_fa_file,
	                                                              aln_bam_file)
	    failed = False
	    for cmd in (index_cmd, aln_cmd):
		try:
		    if self.debug:
			print cmd
		    subprocess.call(cmd, shell=True)
		except:
		    sys.stderr.write('Failed to run:%s\n' % cmd)
		    failed = True
	    if not failed:
		return pysam.AlignmentFile(aln_bam_file)
	    else:
		return None

	query_target = {}
	new_multi_aligns = {}
	align_jobs = []
	query_spans = {}
	for align, block_matches, clipped_ends in partial_aligns:
	    transcript = None
	    target_seq = None
	    for pos, clipped_seq in clipped_ends.iteritems():
		if len(clipped_seq) <= min_len:
		    print '%s:skip partial because it is too short %s %d' % (align.query,
		                                                             clipped_seq,
		                                                             len(clipped_seq))
		    continue
		
		if pos == 'start':
		    query_spans[align.query] = (1, len(clipped_seq))
		else:
		    query_spans[align.query] = (align.query_len - len(clipped_seq) + 1, align.query_len)
		
		if transcript is None:
		    transcript = get_transcript(block_matches)
		    target_seq = transcript.get_sequence(self.genome_fasta)

		if target_seq is not None:
		    matches = search_by_regex(clipped_seq, target_seq)
		    if matches and len(matches) == 1:
			match = matches[0]
			if target_type == 'genome':
			    target = transcript.chrom
			else:
			    target = align.target
			new_align = find_new_align(align.query, query_spans[align.query], target, match, transcript)
			if new_align is not None:
			    new_multi_aligns[align.query] = [align, new_align]
		    else:
			align_jobs.append((align, 
			                   '%s-%s' % (align.query, pos),
			                   transcript.id,
			                   clipped_seq, 
			                   target_seq))
			query_target[align.query] = transcript.id
				    
	if align_jobs:
	    query_fa_file = '%s/partial_query.fa' % self.working_dir
	    target_fa_file = '%s/partial_target.fa' % self.working_dir
	    partial_bam = run_align(align_jobs, query_fa_file, target_fa_file)
	    if partial_bam is not None:
		original_aligns = dict((job[0].query, job[0]) for job in align_jobs)
		new_aligns = parse_partial_aligns(partial_bam, query_target, query_spans)
		
		for query, align in new_aligns.iteritems():
		    new_multi_aligns[query] = [original_aligns[query]]
		    new_multi_aligns[query].append(align)

		if not self.debug:
		    os.remove(partial_bam.filename)
		    os.remove(query_fa_file)
		    for ff in glob.glob(target_fa_file + '*'):
			os.remove(ff)
		
	return new_multi_aligns
    
    def find_indels(self, align, query_fasta, target_fasta, target_type, min_size=0, min_flanking=0, no_indels=False):
	def extract_flanking_blocks():
	    flanking_blocks = []
	    block_idx = -1
	    for i in range(len(align.sam.cigartuples)):
		op, length = align.sam.cigartuples[i]
		if op == 0:
		    block_idx += 1
		    continue
		if (op == 1 or op == 2 or op == 3) and\
		   i > 0 and i < len(align.sam.cigartuples) - 1 and\
		   align.sam.cigartuples[i - 1][0] == 0 and\
		   align.sam.cigartuples[i + 1][0] == 0:
		    flanking_blocks.append((block_idx, block_idx + 1))
	    return flanking_blocks

	def is_del_possibly_intron(target_breaks):
	    motif = target_fasta.fetch(align.target, target_breaks[0], target_breaks[0] + 2) +\
	          target_fasta.fetch(align.target, target_breaks[1] - 3, target_breaks[1] - 1)
	    if motif.lower() == 'gtag' or motif.lower() == 'ctac':
		return True
	    else:
		return False

	adjs = []
		
	if not re.search('[IDN]\d+', align.cigarstring):
	    return adjs

	for i, j in extract_flanking_blocks():
	    target_breaks = (align.blocks[i][1], align.blocks[j][0])
	    target_gap = align.blocks[j][0] - align.blocks[i][1] - 1
	    if align.strand == '+':
		query_breaks = (align.query_blocks[i][1], align.query_blocks[j][0])
		query_gap = align.query_blocks[j][0] - align.query_blocks[i][1] - 1
	    else:
		query_breaks = (align.query_blocks[j][0], align.query_blocks[i][1])
		query_gap = align.query_blocks[i][1] - align.query_blocks[j][0] - 1
		
	    event_type = None
	    # del
	    if target_gap > 0 and query_gap == 0:
		event_type = 'del'
		
	    elif target_gap == 0 and query_gap > 0:
		event_type = 'ins'

	    elif target_gap > 0 and query_gap > 0:
		event_type = 'indel'
		
	    if event_type is None:
		continue

	    if target_type == 'genome' and event_type == 'del' and is_del_possibly_intron(target_breaks):
		continue

	    if event_type in ('ins', 'del') and max(target_gap, query_gap) < min_size:
		if self.debug:
		    print '%s: filter out %s - size too small %d' % (align.query,
		                                                     event_type,
		                                                     max(target_gap, query_gap))
		continue

	    if event_type in ('ins', 'del') and\
	       (query_breaks[0] < min_flanking or\
	        align.query_len - query_breaks[1] + 1 < min_flanking):
		if self.debug:
		    print '%s: filter out %s - flanking too small %d, %d' % (align.query,
		                                                             event_type,
		                                                             query_breaks[0],
		                                                             align.query_len - query_breaks[1] + 1)
		continue

	    adj = Adjacency(align.query,
	                    (align.target, align.target),
	                    query_breaks,
	                    target_breaks,
	                    rearrangement = event_type,
	                    orients = ('L', 'R')
	                    )

	    if event_type == 'ins':
		ins_seq = query_fasta.fetch(align.query, query_breaks[0], query_breaks[1] - 1)
		if align.strand == '-':
		    ins_seq = reverse_complement(ins_seq)
		adj.ins_seq = ins_seq
		if not self.is_repeat_number_change(adj,
	                                            query_fasta, 
	                                            target_fasta,
	                                            align.strand):
		    if len(ins_seq) >= 10:
			self.is_duplication(adj,
		                            query_fasta.fetch(align.query),
		                            target_fasta,
		                            align.strand)

	    if event_type == 'del':
		self.is_repeat_number_change(adj,
	                                     query_fasta,
	                                     target_fasta,
	                                     align.strand)

	    if no_indels and adj.rearrangement in ('ins', 'del'):
		continue

	    adjs.append(adj)

	return adjs
    
    def is_duplication(self, adj, query_seq, target_fasta, strand, min_size_to_align=20, min_size_for_terminal_snp=20):
	def update_support_span(dup_is_up_in_target):
	    seq_breaks = sorted(adj.seq_breaks)
	    if dup_is_up_in_target:
		if strand == '+':
		    adj.support_span = (seq_breaks[0], seq_breaks[0] + 1)
		else:
		    adj.support_span = (seq_breaks[1] - 1, seq_breaks[1])
	    else:
		if strand == '+':
		    adj.support_span = (seq_breaks[1] - 1, seq_breaks[1])
		else:
		    adj.support_span = (seq_breaks[0], seq_breaks[0] + 1)

	query_name = adj.seq_id
	target_name = '%s:%s-%s' % (adj.targets[0], adj.target_breaks[0], adj.target_breaks[1])
	
	novel_seq = query_seq[adj.seq_breaks[0] : adj.seq_breaks[1] - 1]
	if self.is_homopolymer(novel_seq):
	    return
	if strand == '-':
	    novel_seq = reverse_complement(novel_seq)
	    	
	try:
	    target_downstream_seq = target_fasta.fetch(adj.targets[0], adj.target_breaks[1] - 1, adj.target_breaks[1] - 1 + len(novel_seq))
	    target_upstream_seq = target_fasta.fetch(adj.targets[0], adj.target_breaks[0] - len(novel_seq), adj.target_breaks[0])
	except:
	    return
    
	novel_seqs = [novel_seq]
	if len(novel_seq) >= min_size_for_terminal_snp:
	    novel_seqs.append(novel_seq[1:])
	    novel_seqs.append(novel_seq[:-1])
	for ns in novel_seqs:
	    matches_downstream = search_by_regex(ns, target_downstream_seq)
	    if matches_downstream:
		break
	if not matches_downstream:
	    for ns in novel_seqs:
		matches_upstream = search_by_regex(ns, target_upstream_seq)
		if matches_upstream:
		    break

	# no matches, try alignment
	#if not matches_downstream and not matches_upstream and len(novel_seq) >= min_size_to_align:
	    #matches_downstream = search_by_align(novel_seq, target_downstream_seq, query_name, target_name, self.working_dir, debug=self.debug)
	    #if not matches_downstream:
		#matches_upstream = search_by_align(novel_seq, target_upstream_seq, query_name, target_name, self.working_dir, debug=self.debug)

	if matches_downstream:
	    matches = matches_downstream
	    if matches[0][0] == 1:
		# duplication
		if len(matches) == 1:
		    adj.rearrangement = 'dup'
		    adj.target_breaks = (adj.target_breaks[1] + len(novel_seq) - 1, adj.target_breaks[1])
		    adj.orients = ('L', 'R')
		    update_support_span(False)
		else:
		    # more work needed, novel sequence may actually be more than one copy
		    print 'repeat_expansion?', adj.seq_id, novel_seq
	    else:
		print 'gap between dup %s %s' % (adj.seq_id, novel_seq)
		
	elif matches_upstream:
	    matches = matches_upstream
	    if matches[-1][1] == len(target_upstream_seq):
		if len(matches) == 1:
		    adj.rearrangement = 'dup'
		    adj.target_breaks = (adj.target_breaks[0], adj.target_breaks[0] - len(novel_seq) + 1)
		    adj.orients = ('L', 'R')
		    update_support_span(True)
		else:
		    # more work needed, novel sequence may actually be more than one copy
		    print 'repeat_expansion?', adj.seq_id, novel_seq
	    else:
		print 'gap between dup %s %s' % (adj.seq_id, novel_seq)

    def adjust_for_amino_acid_repeat(self, adj):
	def get_aa_repeat(aa, nt, start=False, end=False):
	    repeat = None
	    same = True
	    if not aa:
		return '', ''

	    elif start:
		repeat = aa[0]
		repeat_nt = nt[:3]
		i = 1
		j = 3
		while same and i < len(aa):
		    if aa[i] == repeat[0] and nt[j:j+3] == repeat_nt[:3]:
			repeat += aa[i]
			repeat_nt += nt[j:j+3]
			same = True
		    else:
			same = False
		    i += 1
		    j += 3
	    elif end:
		repeat = aa[-1]
		repeat_nt = nt[-3:]
		i = len(aa) - 2
		j = len(nt) - 6
		while same and i >= 0:
		    if aa[i] == repeat[0] and nt[j:j+3] == repeat_nt[:3]:
			repeat += aa[i]
			repeat_nt += nt[j:j+3] 
			same = True
		    else:
			same = False
		    i -= 1
		    j -= 3
	    return repeat, repeat_nt.lower()

	def update_adj(genome_breaks, repeat_seq, copy_num_change):
	    adj.genome_breaks = genome_breaks
	    adj.repeat_seq = repeat_seq.upper()
	    adj.copy_num_change = copy_num_change

	def construct_seq(chrom, span, break_start, break_end=None, repeat_size=3, add=False, minus=False, repeat_seq=None, flank=50):
	    """start = begining coord of repeat, end = end coord

	    In some cases an event may come from a genomic contig and lie close to an exon boundary,
	    this may cause the conversion from transcript to genome coordinate to not make sense
	    and the sequence asked to be constructed impossible to do - in that case, None will be reported

	    """
	    if add and repeat_seq is not None:
		try:
		    return self.genome_fasta.fetch(chrom, span[0] - flank - 1, break_start).lower() +\
		           repeat_seq.upper() +\
		           self.genome_fasta.fetch(chrom, break_start, span[1] + flank).lower()
		except:
		    return None
	    elif minus and break_end is not None:
		try:
		    return self.genome_fasta.fetch(chrom, span[0] - flank - 1, break_start - 1).upper() +\
		           self.genome_fasta.fetch(chrom, break_end, span[1] + flank).lower()
		except:
		    return None

	def identify_aa_repeat():
	    def search_for_repeat(aa, pos):
		for start in (start - 1, start + 3):
		    before = aa[pos - 1]
		    after = aa[pos]

	    aa = adj.transcripts[0].translate(self.genome_fasta)
	    cds = adj.transcripts[0].get_sequence(self.genome_fasta, cds_only=True)

	    mid_genome_break = ((adj.genome_breaks[0] + 1) + (adj.genome_breaks[1] - 1)) / 2
	    cds_pos = adj.transcripts[0].genome_coord_to_txt_coord(mid_genome_break, cds=True)
	    repeat_spans = []
	    if cds_pos is not None:
		aa_pos = adj.transcripts[0].txt_coord_to_aa_coord(cds_pos)
		if adj.genome_breaks[1] - adj.genome_breaks[0] < 4:
		    aa_pos -= 1
		cds_pos_converted = adj.transcripts[0].aa_coord_to_txt_coord(aa_pos)[1]
		repeat_down, repeat_down_nt = get_aa_repeat(aa[aa_pos:], cds[cds_pos_converted:], start=True)
		repeat_up, repeat_up_nt = get_aa_repeat(aa[:aa_pos], cds[:cds_pos_converted], end=True)

		# identify amino acid repeat start and end
		repeat_start = None
		repeat_end = None
		if len(repeat_down) > 1 or len(repeat_up) > 1:
		    shift_result = None
		    if len(repeat_down) >= 1 and len(repeat_up) >= 1:
			if repeat_down[0] == repeat_up[0] and repeat_down_nt[:3] == repeat_up_nt[:3]:
			    repeat_end = aa_pos + 1 + (len(repeat_down) - 1)
			    repeat_start = aa_pos - (len(repeat_up) - 1)
			    repeat_spans.append((repeat_start, repeat_end))

			else:
			    if len(repeat_down) > 1:
				repeat_end = aa_pos + 1 + (len(repeat_down) - 1)
				repeat_start = aa_pos + 1
				repeat_spans.append((repeat_start, repeat_end))

			    if len(repeat_up) > 1:
				repeat_start = aa_pos - (len(repeat_up) - 1)
				repeat_end = aa_pos
				repeat_spans.append((repeat_start, repeat_end))

		    elif len(repeat_down) > 1:
			repeat_end = aa_pos + 1 + (len(repeat_down) - 1)
			repeat_start = aa_pos + 1
			repeat_spans.append((repeat_start, repeat_end))

		    else:
			repeat_start = aa_pos - (len(repeat_up) - 1)
			repeat_end = aa_pos
			repeat_spans.append((repeat_start, repeat_end))

		    return repeat_spans

	    return []

	def shift_coord(repeat_start, repeat_end):
	    repeat_track_genome = (adj.transcripts[0].aa_coord_to_genome_coord(repeat_start),
	                           adj.transcripts[0].aa_coord_to_genome_coord(repeat_end))
	    repeat_track_sorted = sorted([repeat_track_genome[0][0],
	                                  repeat_track_genome[0][1],
	                                  repeat_track_genome[1][0],
	                                  repeat_track_genome[1][1]])
	    repeat_start, repeat_end = repeat_track_sorted[0], repeat_track_sorted[-1]
	    repeat_start_genome_sorted = sorted(repeat_track_genome[0])
	    repeat_seq = self.genome_fasta.fetch(adj.chroms[0],
	                                         repeat_start_genome_sorted[0] - 1,
	                                         repeat_start_genome_sorted[1])
	    copy_num_ref = (repeat_end - repeat_start + 1) / len(repeat_seq)
	    copy_num_change = adj.size / len(repeat_seq)
	    if adj.event == 'repeat_expansion':
		old_seq = construct_seq(adj.chroms[0],
		                        (repeat_start, repeat_end),
		                        break_start = adj.genome_breaks[0],
		                        repeat_seq = adj.repeat_seq,
		                        add = True)
		if adj.transcripts[0].strand == '+':
		    break_start = repeat_end
		else:
		    break_start = repeat_start - 1
		new_seq = construct_seq(adj.chroms[0],
		                        (repeat_start, repeat_end),
		                        break_start = break_start,
		                        repeat_seq = repeat_seq * copy_num_change,
		                        add = True)
		if old_seq is not None and new_seq is not None and old_seq.lower() == new_seq.lower():
		    return (repeat_start, repeat_end, repeat_seq, (break_start, break_start + 1), (copy_num_ref, copy_num_ref + copy_num_change))

	    elif adj.event == 'repeat_reduction':
		old_seq = construct_seq(adj.chroms[0],
	                                (repeat_start, repeat_end),
	                                break_start = adj.genome_breaks[0] + 1,
	                                break_end = adj.genome_breaks[1] - 1,
	                                minus = True)
		if adj.transcripts[0].strand == '+':
		    break_start = repeat_end - adj.size + 1
		    break_end = repeat_end
		else:
		    break_start = repeat_start
		    break_end = repeat_start + adj.size - 1
		new_seq = construct_seq(adj.chroms[0],
	                                (repeat_start, repeat_end),
		                        break_start = break_start,
		                        break_end = break_end,
	                                minus = True)
		if old_seq is not None and new_seq is not None and old_seq.lower() == new_seq.lower():
		    return (repeat_start, repeat_end, repeat_seq, (break_start, break_end), (copy_num_ref, copy_num_ref - copy_num_change))
	    return None

	aa_repeat_spans = identify_aa_repeat()
	for aa_repeat_span in aa_repeat_spans:
	    shift_result = shift_coord(aa_repeat_span[0], aa_repeat_span[1])
	    if shift_result is not None:
		repeat_start, repeat_end, repeat_seq, breaks, copy_num_change = shift_result
		if adj.event == 'repeat_expansion':
		    update_adj(breaks, repeat_seq, copy_num_change)
		elif adj.event == 'repeat_reduction':
		    update_adj((breaks[0] - 1, breaks[1] + 1), repeat_seq, copy_num_change)
		break

    def is_repeat_number_change(self, adj, query_fasta, target_fasta, strand, max_size=30):
	def chop_seq(seq, size):
	    """chop seq into equal size sub-strings"""
	    return map(''.join, zip(*[iter(seq)] * size))

	def find_repeat(seq):
	    if len(seq) % 3 == 0:
		uniq_subseqs = Set(chop_seq(seq, 3))
		if len(uniq_subseqs) == 1:
		    return list(uniq_subseqs)[0]

	    if len(seq) % 2 == 0:
		uniq_subseqs = Set(chop_seq(seq, 2))
		if len(uniq_subseqs) == 1:
		    return list(uniq_subseqs)[0]
	    return None
	
	def get_copy_num(target, start_pos, repeat, direction):
	    pos = start_pos
	    copy_num = 0
	    matched = True
	    while matched:
		if direction == '+':
		    start = pos - 1
		    end = pos - 1 + len(repeat)
		else:
		    start = pos - len(repeat)
		    end = pos
    
		next_copy = target_fasta.fetch(target, start, end)
		if next_copy.lower() == repeat.lower():
		    copy_num += 1
		    if direction == '+':
			pos = end + 1
		    else:
			pos = start
		else:
		    matched = False

	    return copy_num
	
	if adj.rearrangement == 'ins':
	    seq_breaks_sorted = sorted(adj.seq_breaks)
	    changed_seq = query_fasta.fetch(adj.seq_id, seq_breaks_sorted[0], seq_breaks_sorted[1] - 1)
	    if strand == '-':
		changed_seq = reverse_complement(changed_seq)
	elif adj.rearrangement == 'del':
	    target_breaks_sorted = sorted(adj.target_breaks)
	    # don't need to reverse-complement because target sequence is used
	    changed_seq = target_fasta.fetch(adj.targets[0], target_breaks_sorted[0], target_breaks_sorted[1] - 1)	    
	
	if len(changed_seq) <= max_size and (len(changed_seq) % 2 == 0 or len(changed_seq) % 3 == 0) and not self.is_homopolymer(changed_seq):
	    repeat = find_repeat(changed_seq)
	    if repeat is not None:
		changed_copy_num = len(changed_seq) / len(repeat)
		target_breaks_sorted = sorted(adj.target_breaks)
		copy_num_upstream = get_copy_num(adj.targets[0], target_breaks_sorted[0], repeat, '-')
		copy_num_downstream = get_copy_num(adj.targets[0], target_breaks_sorted[1], repeat, '+')

		if copy_num_upstream > 0 or copy_num_downstream > 0:
		    original_copy_num = copy_num_upstream +  copy_num_downstream
    
		    # just a plain duplication or deletion
		    if original_copy_num == 1 and changed_copy_num == 1:
			return False
 
		    if adj.rearrangement == 'ins':
			adj.event = 'repeat_expansion'
			adj.copy_num_change = (original_copy_num, original_copy_num + changed_copy_num)
		    elif adj.rearrangement == 'del':
			original_copy_num += changed_copy_num
			adj.event = 'repeat_reduction'
			adj.copy_num_change = (original_copy_num, original_copy_num - changed_copy_num)
		    adj.repeat_seq = repeat.upper()

		    expansion_downstream = len(repeat) * copy_num_downstream
		    expansion_upstream = len(repeat) * copy_num_upstream
		    seq_breaks_sorted = sorted(adj.seq_breaks)
		    if strand == '+':
			adj.support_span = (max(1, seq_breaks_sorted[0] - expansion_upstream),
			                    min(query_fasta.fetch(adj.seq_id), seq_breaks_sorted[1] + expansion_downstream))
		    else:
			adj.support_span = (max(1, seq_breaks_sorted[0] - expansion_downstream),
			                    min(query_fasta.fetch(adj.seq_id), seq_breaks_sorted[1] + expansion_upstream))

		    return True

	return False

    def is_fusion(self, adj, align_strands, target_type):
	if adj.transcripts[0].gene != adj.transcripts[1].gene:
	    adj.event = 'fusion'

	    # if there is homol sequence, and one or both
	    # of both breakpoints is not exon_bound, can we adjust
	    # breakpoint to make it exon_bound
	    if target_type == 'transcripts' and\
	       adj.homol_seq and\
	       (not adj.exon_bounds[0] or not adj.exon_bounds[1]):
		adj.adjust_transcript_breaks(align_strands)

	    self.is_sense_fusion(adj, align_strands, target_type)

	    if adj.sense_fusion and adj.chroms[0] == adj.chroms[1]:
		self.is_fusion_read_through(adj)

	    return True
	return False
	    
    def is_sense_fusion(self, adj, align_strands, target_type):
	def reverse_tuple(tup):
	    return tuple(list(tup)[::-1])

	def set_orients(reverse=False):
	    if not reverse:
		adj.upstream_transcript = adj.transcripts[0]
		adj.downstream_transcript = adj.transcripts[1]
		adj.exons_oriented = adj.exons
		adj.exon_bounds_oriented = adj.exon_bounds
	    else:
		adj.upstream_transcript = adj.transcripts[1]
		adj.downstream_transcript = adj.transcripts[0]
		adj.exons_oriented = reverse_tuple(adj.exons)
		adj.exon_bounds_oriented = reverse_tuple(adj.exon_bounds)

	sense = None
	if (adj.event == 'fusion' or adj.event == 'read_through'):
	    if target_type == 'genome':
		sense = False
		if ((adj.transcripts[0].strand == '+' and adj.orients[0] == 'L') or\
		    (adj.transcripts[0].strand == '-' and adj.orients[0] == 'R'))\
		   and\
		   ((adj.transcripts[1].strand == '+' and adj.orients[1] == 'R') or\
		    (adj.transcripts[1].strand == '-' and adj.orients[1] == 'L')):
		    sense = True
		    set_orients(reverse=False)

		elif ((adj.transcripts[1].strand == '+' and adj.orients[1] == 'L') or\
		      (adj.transcripts[1].strand == '-' and adj.orients[1] == 'R'))\
		     and\
		     ((adj.transcripts[0].strand == '+' and adj.orients[0] == 'R') or\
		      (adj.transcripts[0].strand == '-' and adj.orients[0] == 'L')):
		    sense = True
		    set_orients(reverse=True)

	    elif target_type == 'transcripts':
		sense = False
		if align_strands[0] == align_strands[1]:
		    sense = True
		    if align_strands[0] == '+':
			set_orients(reverse=False)
		    else:
			set_orients(reverse=True)

	    # does not allow downstream breakpoint to be at utr when upstream isn't
	    within_utrs = (adj.transcripts[0].within_utr(adj.genome_breaks[0]),
	                   adj.transcripts[1].within_utr(adj.genome_breaks[1]))

	adj.sense_fusion = sense

    def is_fusion_read_through(self, adj):
	"""Check if fusion event is read_through

	   Criteria:
	   1. sense fusion on same chromosome
	   2a. the 2 genes overlap
	   2b. there is no coding gene on the same strand between the 2 genes
	"""
	if adj.chroms[0] != adj.chroms[1] or adj.transcripts[0].strand != adj.transcripts[1].strand:
	    return

	# make sure upstream and downstream gene make sense
	if (adj.upstream_transcript.strand == '+' and adj.upstream_transcript.exons[0][0] > adj.downstream_transcript.exons[0][0]) or\
	   (adj.upstream_transcript.strand == '-' and adj.downstream_transcript.exons[-1][1] > adj.upstream_transcript.exons[-1][1]):
	    return

	# find interval between genes
	interval = adj.transcripts[0].exons[-1][1] + 1, adj.transcripts[1].exons[0][0] - 1
	if adj.transcripts[0].strand == '+':
	    if adj.transcripts[0].exons[0][0] < adj.transcripts[1].exons[0][0]:
		transcripts = adj.transcripts
	    else:
		transcripts = adj.transcripts[::-1]
	    interval = transcripts[0].exons[-1][1] + 1, transcripts[1].exons[0][0] - 1
	if adj.transcripts[0].strand == '-':
	    if adj.transcripts[0].exons[-1][1] > adj.transcripts[1].exons[-1][1]:
		transcripts = adj.transcripts
	    else:
		transcripts = adj.transcripts[::-1]
	    interval = transcripts[1].exons[-1][1] + 1, transcripts[0].exons[0][0] - 1
	
	# if gene1 and gene2 overlaps
	if not interval[0] < interval[1]:
	    adj.event = 'read_through'

	# check if there is any coding gene on same strand within interval
	else:
	    genes_in_between = Set()
	    gene_previous = None
	    for feature in self.annot.fetch(adj.chroms[0], interval[0], interval[1]):
		if not 'exon' in str(feature):
		    continue

		match = re.search(r'transcript_id "(\S+)"', str(feature))
		if match:
		    transcript_id = match.group(1)
		    transcript = self.transcripts_dict[transcript_id]
		    if gene_previous is not None:
			if transcript.gene == gene_previous:
			    continue
		    else:
			gene_previous = transcript.gene
		    if transcript.gene != adj.transcripts[0].gene and\
		       transcript.gene != adj.transcripts[1].gene and\
		       transcript.strand == adj.transcripts[0].strand and\
		       transcript.is_coding() and\
		       transcript.exons[0][0] > interval[0] and\
		       transcript.exons[-1][1] < interval[1]:
			genes_in_between.add(transcript.gene)
			break

	    if not genes_in_between:
		adj.event = 'read_through'
    
    def is_read_through_from_single_align(self, block_matches, align):
	def check_matches(matches, index, before=False, after=False):
	    """checks every match before or after an index is None"""
	    if before:
		unmatched = True
		for i in range(0, index):
		    if matches[i] is not None:
			unmatched = False
			break
		return unmatched
	    elif after:
		unmatched = True
		for i in range(index + 1, len(matches)):
		    if matches[i] is not None:
			unmatched = False
			break
		return unmatched
	    else:
		return None

	genes = Set()
	for tid in block_matches.keys():
	    if self.transcripts_dict.has_key(tid):
		genes.add(self.transcripts_dict[tid].gene)
		
	if len(genes) < 2:
	    return None

	upstream_matches = {}
	downstream_matches = {}
	for tid, matches in block_matches.iteritems():
	    for i in range(len(matches) - 1):
		if matches[i] is None and\
		   matches[i + 1] is not None and\
		   matches[i + 1][0][1][0] == '=' and\
		   check_matches(matches, i, before=True):
		        
		    if not downstream_matches.has_key(i + 1):
			downstream_matches[i + 1] = []
		    downstream_matches[i + 1].append((tid, matches[i + 1][0][0], len([m for m in matches if m is not None])))
		elif matches[i + 1] is None and\
		     matches[i] is not None and \
		     matches[i][-1][1][1] == '=' and\
		     check_matches(matches, i + 1, after=True):
		    if not upstream_matches.has_key(i):
			upstream_matches[i] = []
		    upstream_matches[i].append((tid, matches[i][-1][0], len([m for m in matches if m is not None])))
		    
	if len(upstream_matches.keys()) == 1 and len(downstream_matches.keys()) == 1:
	    upstream_block_idx = upstream_matches.keys()[0]
	    downstream_block_idx = downstream_matches.keys()[0]
	    
	    # allows novel exon in the junction
	    if downstream_block_idx > upstream_block_idx:
		upstream_matches_sorted = sorted(upstream_matches[upstream_block_idx], key = itemgetter(2), reverse=True)
		downstream_matches_sorted = sorted(downstream_matches[downstream_block_idx], key = itemgetter(2), reverse=True)
		upstream_transcript = self.transcripts_dict[upstream_matches_sorted[0][0]]
		downstream_transcript = self.transcripts_dict[downstream_matches_sorted[0][0]]
		
		upstream_gene = upstream_transcript.gene
		downstream_gene = downstream_transcript.gene
		
		upstream_exon = upstream_transcript.exon_num(upstream_matches_sorted[0][1])
		downstream_exon = downstream_transcript.exon_num(downstream_matches_sorted[0][1])
		exons = [upstream_exon, downstream_exon]
		seq_breaks = [align.query_blocks[upstream_block_idx][1], align.query_blocks[downstream_block_idx][0]]
		target_breaks = [align.blocks[upstream_block_idx][1], align.blocks[downstream_block_idx][0]]
		orients = ['L', 'R']

		adj = Adjacency(align.query,
		                (align.target, align.target),
		                seq_breaks,
		                target_breaks,
		                event = 'read_through',
		                #transcripts = transcripts,
		                exons = exons,
		                exon_bounds = (True, True),
		                chroms = (align.target, align.target),
		                genome_breaks = target_breaks,
		                orients = orients,
		                sense_fusion = True,
		                )
		
		adj.transcripts = self.map_transcripts_to_adj(adj, (align, align), block_matches=block_matches,
		                                              genes = (upstream_gene, downstream_gene))
		if adj.transcripts[0].strand == '+':
		    adj.upstream_transcript = adj.transcripts[0]
		    adj.downstream_transcript = adj.transcripts[1]
		    adj.exons_oriented = adj.exons
		    adj.exon_bounds_oriented = adj.exon_bounds
		else:
		    adj.upstream_transcript = adj.transcripts[1]
		    adj.downstream_transcript = adj.transcripts[0]
		    adj.exons_oriented = tuple(reversed(list(adj.exons)))
		    adj.exon_bounds_oriented = tuple(reversed(list(adj.exon_bounds)))
		adj.update_transcript_breaks()
		
		return adj
	return None
    
    @classmethod
    def set_frame(cls, adjs, query_fasta, genome_fasta):
	for adj in adjs:
	    # already set in repeat_expansion/reduction cases
	    if adj.in_frame:
		continue

	    seq_id = adj.seq_id.split(',')[0]
	    query_seq = query_fasta.fetch(seq_id)
	    seq_breaks = sorted(map(int, adj.seq_breaks.split(',')[0].split('-')))
	    genome_breaks = sorted(map(int, list(adj.genome_breaks)))

	    transcripts = adj.transcripts
	    if adj.upstream_transcript and adj.downstream_transcript:
		transcripts = (adj.upstream_transcript, adj.downstream_transcript)
		
	    if transcripts:
		in_frame = check_frame(query_seq,
		                       seq_breaks,
		                       genome_breaks, 
		                       genome_fasta,
		                       transcripts,
		                       adj.event,
		                       adj.size)
		if type(in_frame) is tuple:
		    adj.in_frame = True
		else:
		    adj.in_frame = in_frame
