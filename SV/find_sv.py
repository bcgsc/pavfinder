from optparse import OptionParser
from itertools import groupby
from sets import Set
import sys
import os
import re
import pysam
import gapped_align
import split_align
from shared import gmap, bwa_mem
from variant import Adjacency, Variant
from shared.annotate import overlap_pe, parallel_parse_overlaps, annotate_rna_event, annotate_gene_fusion, update_features, get_acen_coords
from shared.read_support import scan_all, fetch_support
from shared.alignment import reverse_complement, target_non_canonical
from vcf import VCF
from pybedtools import BedTool

class SVFinder:
    
    name = 'PAVFinder'
    version = '0.0.1'
    
    def __init__(self, bam_file, aligner, contig_fasta, genome_fasta, out_dir,
                 genome=None, index_dir=None, num_procs=0, is_transcriptome=None,
                 annot_file=None, skip_simple_repeats=False, cytobands_file=None, acen_buffer=0, debug=False):
	self.bam = pysam.Samfile(bam_file, 'rb')
	self.aligner = aligner
	self.contig_fasta_file = contig_fasta
	self.contig_fasta = pysam.Fastafile(contig_fasta)
	self.ref_fasta = pysam.Fastafile(genome_fasta)
	self.genome = genome
	self.index_dir = index_dir
	self.num_procs = num_procs
	self.is_transcriptome = is_transcriptome
	self.out_dir = out_dir
	self.annot_file = annot_file
	self.adjs = []
	self.skip_simple_repeats = skip_simple_repeats
	self.cytobands_file = cytobands_file
	self.acen_buffer = acen_buffer
	self.debug = debug
	
	self.avg_tlen = None
	self.avg_tlen_normal = None
	
    def find_adjs(self, min_ctg_cov, max_size=None, min_size=None, ins_as_ins=False, skip_acen=False, check_alt_paths=False, min_ctg_size=0, bad_coords=None, skip_contigs_file=None):
	"""Main method to go through the BAM file, extract split and gapped alignments, and calls
	the respective modules to identify adjs"""
	def find_events_in_single_align(align):
	    """Implement as sub-function so that small-scale events can be found on split alignments too"""
	    adjs = gapped_align.find_adjs(align, contig_seq, self.is_transcriptome, ins_as_ins=ins_as_ins)
		
	    repeats = Set()
	    for i in range(len(adjs)):
		adj = adjs[i]
		
		if self.skip_simple_repeats and self.break_region_has_low_complexity(adj.chroms[0], adj.breaks):
		    repeats.add(i)
		    if self.debug:
			sys.stdout.write("remove contig %s %s potential simple-repeat %s:%s-%s\n" % (adj.contigs[0], 
			                                                                             adj.rearrangement, 
			                                                                             adj.chroms[0], 
			                                                                             adj.breaks[0], 
			                                                                             adj.breaks[1]))
		    continue
		    
		new_contig_breaks = self.expand_contig_breaks(adj.chroms[0], adj.breaks, contig, adj.contig_breaks[0], adj.rearrangement, self.debug)
		if new_contig_breaks is not None:
		    adj.contig_breaks[0] = new_contig_breaks
		    			    
	    if repeats:
		for i in sorted(repeats, reverse=True):
		    del adjs[i]
		
	    return adjs	
	
	def is_align_in_acen(align, acen):
	    """Checks to see if alignment overlaps with acentromeric coordinates
	    Args:
	        align: alignment (Alignment)
		acen: acentromeric coordinates parsed from UCSC cytobands file (Dictionary) {chrom:(start, end), (start, end)}
	    Returns True if overlapped
	    """
	    s1, e1 = align.tstart, align.tend
	    if acen.has_key(align.target):
		for (start, end) in acen[align.target]:
		    s2, e2 = int(start) - self.acen_buffer, int(end) + self.acen_buffer
		    if s1 <= e2 and s2 <= e1:
			return True
		    
	    return False
	
	def create_set(list_file):
	    subset = Set()
	    for line in open(list_file, 'r'):
		subset.add(line.strip('\n'))
	    return subset
		    
	acen_coords = None
	if skip_acen:
	    acen_coords = get_acen_coords(self.cytobands_file)
	    
	skip_contigs = None
	if skip_contigs_file and os.path.exists(skip_contigs_file):
	    skip_contigs = create_set(skip_contigs_file)
	
	all_adjs = []
	for contig, group in groupby(self.bam.fetch(until_eof=True), lambda x: x.qname):
	    print 'contig', contig
	    alns = list(group)
	    contig_seq = self.contig_fasta.fetch(contig)
	    
	    if len(contig_seq) < min_ctg_size:
		if self.debug:
		    sys.stdout.write('%s(%d bp) less than min contig size %d bp\n' % (contig, len(contig_seq), min_ctg_size))
		continue
	    
	    if skip_contigs and contig in skip_contigs:
		if self.debug:
		    sys.stdout.write('%s skipped\n' % contig)
		continue
	    
	    if len(alns) > 1:
		chimeric_aligns, dubious = split_align.find_chimera(alns, 
		                                                    self.aligner, 
		                                                    self.bam, 
		                                                    min_coverage=min_ctg_cov, 
		                                                    check_alt_paths=check_alt_paths, 
		                                                    debug=self.debug)		
		if chimeric_aligns:
		    if acen_coords:
			skip = False
			for align in chimeric_aligns:
			    if acen_coords and is_align_in_acen(align, acen_coords):
				if self.debug:
				    sys.stdout.write('skip contig %s because alignment is in centromere %s:%d-%d\n' % (contig,
					                                                                               align.target,
					                                                                               align.tstart,
					                                                                               align.tend
					                                                                               ))
				skip = True
				break
			if skip:
			    continue
		
		    adjs = split_align.find_adjs(chimeric_aligns, self.aligner, contig_seq, dubious=dubious, debug=self.debug)
		    
		    bad = Set()
		    for i in range(len(adjs)):
			adj = adjs[i]
			#check if homol is simple repeat
			if adj.homol_seq and adj.homol_seq[0] != '-' and self.is_homol_low_complexity(adj):
			    if self.debug:
				sys.stderr.write("homol_seq is simple-repeat %s:%s\n" % (adj.contigs[0], adj.homol_seq[0]))
			    bad.add(i)
			
			# check if event is simple repeat expansions
			if self.skip_simple_repeats and self.is_novel_sequence_repeat(adj):
			    if self.debug:
				sys.stderr.write("novel_seq is simple-repeat %s:%s\n" % (adj.contigs[0], adj.novel_seq))
			    bad.add(i)
			    
			if i > 0:
			    if adjs[i].chroms == adjs[i - 1].chroms and\
			       adjs[i].breaks == adjs[i - 1].breaks and\
			       adjs[i].orients == adjs[i].orients and\
			       adjs[i].contig_breaks != adjs[i - 1].contig_breaks:
				if self.debug:
				    sys.stderr.write("%s has 2 contig_breaks for same event\n" % adj.contigs[0])
				bad.add(i - 1)
				bad.add(i)
		    
		    if bad:
			for i in sorted(bad, reverse=True):
			    del adjs[i]
			    
		    all_adjs.extend(adjs)
		   
		    # capture small-scale events within each chimeric alignment
		    for align in chimeric_aligns:
			all_adjs.extend(find_events_in_single_align(align))
			
	    best_align = gapped_align.find_single_unique(alns, self.aligner, self.bam, debug=self.debug)
	    if best_align:
		all_adjs.extend(find_events_in_single_align(best_align))
		    
	merged_adjs = Adjacency.merge(all_adjs)
	
	# screen out adjacencies that overlap segdups
	if bad_coords is not None and os.path.exists(bad_coords):
	    self.screen_by_coordinate(merged_adjs, bad_coords)
	
	if max_size is not None:
	    selected = []
	    for adj in merged_adjs:
		size = adj.get_size()
		if size is not None and size != 'NA' and size <= max_size:
		    selected.append(adj)
	    return selected
	elif min_size is not None:
	    selected = []
	    for adj in merged_adjs:
		size = adj.get_size()
		if size is None or size == 'NA' or size >= min_size:
		    selected.append(adj)
	    return selected
	else:
	    return merged_adjs
    
    def create_variants(self, adjs):
	self.variants = []
	# special cases for imprecise insertions
	ins_variants, ins_adjs = Adjacency.extract_imprecise_ins([adj for adj in adjs if adj.align_types[0] == 'split' and adj.rearrangement != 'inv'], debug=self.debug)
	if ins_variants:
	    self.variants.extend(ins_variants)
	    		
	# collect adjs that don't require grouping
	adjs_no_grouping = [adj for adj in adjs if adj.rearrangement != 'inv' and adj.rearrangement != 'trl' and not adj.id in ins_adjs]
	
	# handle inversions
	invs = [adj for adj in adjs if adj.rearrangement == 'inv']
	if invs:
	    self.variants.extend(Adjacency.group_inversions(invs))
	    	    
	# convert translocations to insertions
	trls = [adj for adj in adjs if adj.rearrangement == 'trl']
	trls_remained = None
	if trls:
	    ins_variants, trls_remained = Adjacency.extract_interchrom_ins(trls)
	    self.variants.extend(ins_variants)
	    	    
	# group reciprocal transcloations
	if trls_remained:
	    self.variants.extend(Adjacency.group_trls(trls_remained))
	    	    
	for adj in adjs_no_grouping:
	    if not adj.dubious:
		self.variants.append(Variant(adj.rearrangement.upper(), [adj]))
		
		
    def is_novel_sequence_repeat(self, adj, min_len=3):
	"""Check to see if novel sequence is part of a tandem repeat"""
	contig_seq = self.contig_fasta.fetch(adj.contigs[0]).upper()	
	contig_breaks = (adj.contig_breaks[0][0], adj.contig_breaks[0][1] - 2)

	len_seq = 0
	if adj.novel_seq and len(adj.novel_seq) > min_len:
	    len_seq = len(adj.novel_seq)
	copies = []
	
	if len_seq > 0:
	    for i in xrange(len(contig_seq) - len(adj.novel_seq) + 1):
		if i == contig_breaks[0]:
		    continue
		
		if contig_seq[i:i + len_seq] == adj.novel_seq:
		    # check if it overlaps (or contiguous) with contig_breaks		
		    if min(i + len(adj.novel_seq) - 1, contig_breaks[1]) - max(i, contig_breaks[0]) >= -1:
			copies.append(i)
		    
	return len(copies) > 0
    
    def break_region_has_low_complexity(self, chrom, breaks, min_units=4, buf=1):
	r = re.compile(r"(.+?)\1+")
	window = 100

	for i in range(breaks[0] - window, breaks[1] + 1):
	    try:
		seq = self.ref_fasta.fetch(chrom, max(0,i), i + 100)
	    except:
		print "can't extract reference sequence for complexity checking %s:%s-%s" % (chrom, i, i + 100)
		continue
	    
	    repeats = r.findall(seq)
	    
	    if repeats and repeats[0].upper() != 'N':
		m = re.search('^(%s){1,}' % repeats[0], seq)
		if m is not None:
		    repeat_start = i + 1
		    repeat_end = i + len(m.group(0))
		    num_units = len(m.group(0)) / len(repeats[0])
		    if num_units >= min_units:
			if breaks[0] != breaks[1] and\
			   breaks[0] + 1 >= repeat_start - buf and breaks[1] - 1 <= repeat_end + buf:
			    if self.debug:
				sys.stdout.write('%s %s:%s-%s in low-complexity region %s:%s-%s %sx%d\n' % ('del',
				                                                                            chrom,
				                                                                            breaks[0],
				                                                                            breaks[1],
				                                                                            chrom,
				                                                                            repeat_start,
				                                                                            repeat_end,
				                                                                            repeats[0].upper(),
				                                                                            num_units,
				                                                                            ))
			    
			    return True
			elif breaks[0] == breaks[1] and\
			     breaks[0] >= repeat_start - buf and breaks[0] + 1 <= repeat_end + buf:
			    if self.debug:
				sys.stdout.write('%s %s:%s-%s in low-complexity region %s:%s-%s %sx%d\n' % ('ins',
				                                                                            chrom,
				                                                                            breaks[0],
				                                                                            breaks[1],
				                                                                            chrom,
				                                                                            repeat_start,
				                                                                            repeat_end,
				                                                                            repeats[0].upper(),
				                                                                            num_units,
				                                                                            ))
			    
			    
			    return True
			
	return False
	    
    def extract_repeat(self, seq):
	repeat = {'start':None, 'end':None}
	if len(seq) == 1:
	    repeat['start'] = seq
	    repeat['end'] = seq
	else:
	    re_start = re.compile(r"^(.+?)\1+")
	    re_end = re.compile(r"(.+?)\1+$")

	    repeats = re_start.findall(seq)
	    if repeats:
		repeat['start'] = repeats[0]
	    repeats = re_end.findall(seq)
	    if repeats:
		repeat['end'] = repeats[0]
	
	return repeat
	
    def expand_contig_breaks(self, chrom, breaks, contig, contig_breaks, event, debug=False):
	contig_breaks_sorted = sorted(contig_breaks)
	pos_strand = True if contig_breaks[0] < contig_breaks[1] else False
	contig_seq = self.contig_fasta.fetch(contig)
	contig_breaks_expanded = [contig_breaks[0], contig_breaks[1]]
	
	seq = None
	if event == 'del':
	    seq = self.ref_fasta.fetch(chrom, breaks[0], breaks[1] - 1)
	    if not pos_strand:
		seq = reverse_complement(seq)
		
	elif event == 'ins':
	    seq = contig_seq[contig_breaks_sorted[0] : contig_breaks_sorted[1] - 1]
	        
	if seq is None:
	    print contig, event, 'cannot find seq', seq
	    return None
	
	repeat = self.extract_repeat(seq)
	
	if repeat['end'] is not None:
	    seq = repeat['end']
	    size = len(seq)
	    # downstream
	    start = contig_breaks_sorted[1] - 1
	    expand = 0
	    while start + size <= len(contig_seq):
		next_seq = contig_seq[start : start + size]
		if next_seq.upper() != seq.upper():
		    break
		else:
		    expand += size
		    start += size
	    if pos_strand:
		contig_breaks_expanded[1] += expand
	    else:
		contig_breaks_expanded[0] += expand
	    
	# upstream
	if repeat['start'] is not None:
	    seq = repeat['start']
	    size = len(seq)
	    start = contig_breaks_sorted[0]
	    expand = 0
	    while start - size >= 0:
		next_seq = contig_seq[start - size : start]
		if next_seq.upper() != seq.upper():
		    break
		else:
		    expand -= size
		    start -= size
	    if pos_strand:
		contig_breaks_expanded[0] += expand
	    else:
		contig_breaks_expanded[1] += expand
	    
	if debug and tuple(contig_breaks_expanded) != contig_breaks:
	    sys.stdout.write('contig breaks expanded:%s %s -> %s\n' % (contig, contig_breaks, contig_breaks_expanded))
	    
	return tuple(contig_breaks_expanded)
    	    
    def write_subseq(self, adj, out, name_sep):
	"""Outputs the sub-sequence of a split alignment to output file
	
	Args:
	    adj: Variant object
	    out: Filehandle of output file
	    name_sep: (str) Character used to combined various info into query name
	"""	    
	subseqs = adj.extract_subseqs(self.contig_fasta)
	for i in range(len(subseqs)):
	    out.write('>%s%s%s%s%d\n%s\n' % (adj.contigs[0], name_sep, adj.key(), name_sep, i, subseqs[i]))
	    
    def is_homol_low_complexity(self, adj, min_len=5):
	if adj.homol_seq[0] is not None and len(adj.homol_seq[0]) >= min_len:
	    if self.is_seq_low_complexity(adj.homol_seq[0]):
		return True
	    
	return False
    
    def is_subseq_low_complexity(self, adj):
	subseqs = adj.extract_subseqs(self.contig_fasta)
	
	for subseq in subseqs:
	    if self.is_seq_low_complexity(subseq):
		return True
	    
	return False
    
    def is_seq_low_complexity(self, seq, threshold=0.9):
	bases = ('A', 'G', 'T', 'C')
	
	for base in bases:
	    if float(seq.upper().count(base)) / float(len(seq)) > threshold:
		return True
	    
	    # dimer content
	    for i in xrange(len(bases)-1):
		for j in xrange(i+1, len(bases)):
		    if i == j:
			continue
		    dimer = bases[i] + bases[j]
		    bases_in_dimers = seq.upper().count(dimer) * 2
		    
		    if float(bases_in_dimers) / float(len(seq)) > threshold:
			return True
		
	return False
		
    def write_probe(self, adj, out, name_sep):
	"""Outputs the probe sequence to output file
	
	Args:
	    adj: Variant object
	    out: Filehandle of output file
	    name_sep: (str) Character used to combined various info into query name
	"""
	out.write('>%s%s%s\n%s\n' % (adj.contigs[0], name_sep, adj.key(), adj.probes[0]))
		        
    
    def realign(self, name_sep, use_realigns=False):
	"""Aligns probe and subsequences against genome again
	
	Will generate a Fasta file in the output direcotyr called 'realign.fa' containing all the 
	(subsequences + probe) sequences to be realigned
	Will run the alignment on the command line and generate a BAM file in the output directory 
	called 'realign.bam'
	
	Args:
	    name_sep: (str) Character used to combined various info into query name
	""" 
	prefix = 'realign'
	if not use_realigns:
	    out_file = '%s/%s.fa' % (self.out_dir, prefix)
	    out = open(out_file, 'w')
	    for variant in self.variants:
		for adj in variant.adjs:
		    self.write_probe(adj, out, name_sep)	    
		    #if adj.align_types[0] == 'split':
			#self.write_subseq(adj, out, name_sep)
	    out.close()
	
	# run aligner
	realign_bam_file = '%s/%s.bam' % (self.out_dir, prefix)
	if not use_realigns:
	    if self.aligner == 'gmap':
		return_code = gmap.run(out_file, realign_bam_file, self.genome, self.index_dir, self.num_procs)
	    elif self.aligner == 'bwa_mem':
		return_code = bwa_mem.run(out_file, realign_bam_file, self.genome, self.index_dir, self.num_procs)
	
	return realign_bam_file
    
    
    def screen_realigns(self, use_realigns=False):
	if not self.aligner or not self.genome or not self.index_dir:
	    return None
	
	name_sep = '.'
	realign_bam_file = self.realign(name_sep, use_realigns=use_realigns)
	try:
	    bam = pysam.Samfile(realign_bam_file, 'rb')
	except:
	    sys.exit('Error parsing realignment BAM:%s' % realign_bam_file)
	    
	# creates mapping from query to variant and adjacency
	query_to_variant = {}
	for i in range(len(self.variants)):
	    for j in range(len(self.variants[i].adjs)):
		adj = self.variants[i].adjs[j]
		query = adj.contigs[0] + name_sep + adj.key()
		query_to_variant[query] = (i, j)
		
	failed_variants = Set()
	for key, group in groupby(bam.fetch(until_eof=True), lambda x: name_sep.join(x.qname.split(name_sep)[:2])):
	    alns = list(group)
	    variant_idx = query_to_variant[key][0]
	    variant = self.variants[variant_idx]
	    adj_idx = query_to_variant[key][1]
	    adj = variant.adjs[adj_idx]
	    adj_aligns = adj.aligns[0]
	    
	    indices_to_check = (0, 1)
	    if variant.event == 'INS':
		index = None
		for i in (0, 1):
		    if variant.chrom == adj.chroms[i] and (variant.pos[0] == adj.breaks[i] or variant.pos[1] == adj.breaks[i]):
			index = i
			break
		    
		if index is not None:
		    indices_to_check = (index,)
	    
	    probe_alns = [aln for aln in alns if not aln.qname[-1].isdigit()]
	    if not gapped_align.screen_probe_alns(adj_aligns, probe_alns, adj.align_types[0]):
		if self.debug:
		    sys.stdout.write('probe align completely to one location or not aligned with confidence: %s\n' % key)
		failed_variants.add(variant)
		continue
	    			
	for failed_var in failed_variants:
	    self.variants.remove(failed_var)
		        
    def screen_adjs(self, ins_variants=None, ins_adjs=None):
	"""Screens adjs using realignment results to get rid of non-specific sequences"""
	# screens away potential transcriptome deletions that are actually introns
	if self.is_transcriptome:
	    self.screen_del()
	
	if not self.aligner or not self.genome or not self.index_dir:
	    return None
	
	name_sep = '.'
	realign_bam_file = self.realign(name_sep)
	
	try:
	    bam = pysam.Samfile(realign_bam_file, 'rb')
	except:
	    sys.exit('Error parsing realignment BAM:%s' % realign_bam_file)
        
	adj_dict = dict((self.adjs[i].contigs[0] + name_sep + self.adjs[i].key(), i) for i in range(len(self.adjs)))
	
	# group subsequence alignments by contig name
	# qname of the subsequence alignment is 'contig-index', so split qname and keep the 
	# first element(contig) for grouping
	bad = Set()
	for key, group in groupby(bam.fetch(until_eof=True), lambda x: name_sep.join(x.qname.split(name_sep)[:2])):
	    alns = list(group)
	    adj_idx = adj_dict[key]
	    adj = self.adjs[adj_idx]
	    adj_aligns = adj.aligns[0]
	    	    
	    probe_alns = [aln for aln in alns if not aln.qname[-1].isdigit()]
	    if not gapped_align.screen_probe_alns(adj_aligns, probe_alns):
		if self.debug:
		    sys.stdout.write('probe align completely to one location or not aligned with confidence: %s\n' % key)
		bad.add(adj_idx)
		continue
	    
	    split_alns = [aln for aln in alns if aln.qname[-1].isdigit()]
	    if adj.align_types[0] == 'split' and split_alns:
		if not split_align.screen_subseq_alns(adj_aligns, split_alns, bam, name_sep, debug=self.debug):
		    if self.debug:
			sys.stdout.write('getrid split %s\n' % key)
		    bad.add(adj_idx)
		    		    
	for i in sorted(bad, reverse=True):
	    del self.adjs[i]
	    
    def screen_del(self):
	"""Checks to see if the deletions identified of a transcriptome assembly are actually introns
	Will remove any deletion event if it is potentially an intron
	"""
	bad = Set()
	for i in range(len(self.adjs)):
	    adj = self.adjs[i]
	    if adj.rearrangement != 'del':
		continue
	    	    
	    splice_motif = ''
	    if adj.orients[0] == 'L':
		splice_motif = self.ref_fasta.fetch(adj.chroms[0], adj.breaks[0], adj.breaks[0] + 2)
	    else:
		splice_motif = self.ref_fasta.fetch(adj.chroms[0], adj.breaks[0] - 3, adj.breaks[0] - 1)

	    if adj.orients[1] == 'L':
		splice_motif += self.ref_fasta.fetch(adj.chroms[1], adj.breaks[1], adj.breaks[1] + 2)
	    else:
		splice_motif += self.ref_fasta.fetch(adj.chroms[1], adj.breaks[1] - 3, adj.breaks[1] - 1)
		
	    if splice_motif.upper() == 'GTAG' or reverse_complement(splice_motif).upper() == 'GTAG':
		if self.debug:
		    sys.stdout.write('filtering out deletion %s because it is bounded by splice motif %s\n' % (adj.key(), splice_motif))
		bad.add(i)
		
	for i in sorted(bad, reverse=True):
	    del self.adjs[i]
	    
    def output(self, only_somatic=False, reference_url=None, assembly_url=None, insertion_as_breakends=None):
	variants = [variant for variant in self.variants if not variant.filtered_out]
	if only_somatic:
	    variants = [variant for variant in self.variants if variant.somatic]
	    
	#if variants:
	self.output_variants(variants, 
                             '%s/variants.vcf' % self.out_dir, 
                             reference_url=reference_url, 
                             assembly_url=assembly_url, 
                             insertion_as_breakends=insertion_as_breakends,
	                     )
    
	adjs = []
	for variant in variants:
	    adjs.extend(variant.adjs)
	    
	self.output_adjacencies(adjs, 
                                '%s/adjacencies.tsv' % self.out_dir, 
                                format='tab')
		                         	
    def output_variants(self, variants, out_file, reference_url=None, assembly_url=None, insertion_as_breakends=False):
	"""Output variants in VCF format"""
	records = []
	for variant in variants:
	    output = variant.as_vcf(self.ref_fasta, not self.is_transcriptome, insertion_as_sv=not insertion_as_breakends)
	    
	    if output is not None and output != '':
		records.extend(output.split('\n'))
		
	out = open(out_file, 'w')
	out.write('%s\n' % VCF.header(source='%s_v%s' % (self.name, self.version), reference_url=reference_url, assembly_url=assembly_url))
	# sort by chromosome and pos
	records.sort(key=lambda record: (record.split('\t')[0], record.split('\t')[1]))
	for record in records:
	    out.write('%s\n' % record)
	out.close()
	
    def output_probes(self, adjs, out_file):
	for adj in adjs:
	    for i in range(len(adj.probes)):
		print 'probe', adj.id, i, adj.probes[i]
		
	out = open(out_file, 'w')
	for adj in adjs:
	    if adj.probes[0] != 'NA':
		out.write('>%s %d\n%s\n' % (adj.id, len(adj.probes[0]), adj.probes[0]))
	out.close()
	
    def output_adjacencies(self, adjs, out_file, format):
	"""Output adjacencies in tsv format"""
	fn = None
	args = ()
	if format == 'bedpe':
	    fn = 'as_bedpe'
	elif format == 'tab':
	    fn = 'as_tab'
	    	    
	if not fn is None:
	    out = open(out_file, 'w')
	    if format == 'tab':
		out.write('%s\n' % Adjacency.show_tab_headers())
		
	    for adj in adjs:	    
		output = getattr(adj, fn)(*args)
		try:
		    out.write('%s\n' % output)
		except:
		    sys.stdout.write("can't output adjacency")
		    
	    out.close()
	    	    
    def merge_adjs(self):
	"""Groups contigs bearing same adj together into single adj"""
	self.adjs = Adjacency.merge(self.adjs)
	
    def annotate_adjs(self):
	"""Annotates adjs: genes, transcripts, exons, introns"""
	# output adjs as bedpe
	adjs_bedpe = '%s/adjs.bedpe' % self.out_dir
	self.output_adjs(adjs_bedpe, 'bedpe')
	
	# overlap adj coordinates against annotation file
	overlap_result = '%s/adjs_annot.bed' % self.out_dir
	overlap_pe(adjs_bedpe, self.annot_file, overlap_result)
	
	# update adjs with feature info
	features_overlapped = parallel_parse_overlaps(overlap_result, [adj.key() for adj in self.adjs], self.num_procs)

	for adj in self.adjs:
	    key = adj.key()
	    if not features_overlapped.has_key(key):
		continue
	    
	    update_features(adj, features_overlapped[key])			    
	    annotate_rna_event(adj)
	    if adj.rna_event == 'gene_fusion':
		annotate_gene_fusion(adj)
		
    def find_support(self, tumor_bam=None, min_overlap=None, normal_bam=None, min_overlap_normal=None):
	"""Extracts read support from reads-to-contigs BAM
	
	Assumes reads-to-contigs is NOT multi-mapped 
	(so we add the support of different contigs of the same event together)
	
	Args:
	    bam_file: (str) Path of reads-to-contigs BAM
	"""
	coords = {}
	for variant in self.variants:
	    for adj in variant.adjs:
		for i in range(len(adj.contigs)):
		    contig = adj.contigs[i]
		    span = adj.contig_support_span(i)
				    
		    try:
			coords[contig].append(span)
		    except:
			coords[contig] = [span]
	
	if coords:
	    bam_files = {'tumor':tumor_bam, 'normal':normal_bam}
	    
	    for sample_type in ('tumor', 'normal'):
		bf = bam_files[sample_type]
		if bf is None:
		    continue
		
		if sample_type == 'tumor':
		    avg_tlen = self.avg_tlen
		    overlap_buffer = min_overlap
		    perfect = False
		else:
		    avg_tlen = self.avg_tlen_normal
		    overlap_buffer = min_overlap_normal
		    perfect = False
				    
		avg_tlen = None
		tlens = []
		# if fewer than 10000 adjs, use 'fetch' of Pysam
		if len(coords) < 1000:
		    support, tlens = fetch_support(coords, bf, self.contig_fasta, overlap_buffer=overlap_buffer, perfect=perfect, debug=self.debug)
		# otherwise use multi-process parsing of bam file
		else:
		    support, tlens = scan_all(coords, bf, self.contig_fasta_file, self.num_procs, overlap_buffer=overlap_buffer, perfect=perfect, debug=self.debug)
		    
		print 'total tlens', len(tlens), sample_type
		
		for variant in self.variants:
		    for adj in variant.adjs:
			if sample_type == 'tumor':
			    support_result = adj.support
			    # for sum_support()
			    normal = False
			else:
			    support_result = adj.support_normal
			    normal = True
    
			for i in range(len(adj.contigs)):
			    contig = adj.contigs[i]
			    span = adj.contig_support_span(i)
			    coord = '%s-%s' % (span[0], span[1])
						
			    if support.has_key(contig) and support[contig].has_key(coord):
				support_result['spanning'].append(support[contig][coord][0])
				support_result['flanking'].append(support[contig][coord][1])
				support_result['tiling'].append(support[contig][coord][2])
						    
			adj.sum_support(normal=normal)	
			 
		if tlens:
		    avg_tlen = float(sum(tlens)) / len(tlens)
		    print 'avg tlen', avg_tlen, sample_type
		    	
	if self.debug:
	    for variant in self.variants:
		for adj in variant.adjs:
		    print 'tumor', adj.support, adj.support_total, adj.support_final
		    print 'normal', adj.support_normal, adj.support_total_normal, adj.support_final_normal
				
    def filter_by_read_support(self, adj, min_support, normal=False):
	small = 50
	max_homol_len = 25
	min_flanking = 2
	adj_size = adj.get_size()
	
	if not normal:
	    support_total, avg_tlen = adj.support_total, self.avg_tlen
	else:
	    support_total, avg_tlen = adj.support_total_normal, self.avg_tlen_normal
	    
	support_final = support_total['spanning'] + support_total['flanking']
	filtered_out = False
	
	if support_total['spanning'] < min_support:
	    filtered_out = True
	    if self.debug:
		sys.stdout.write('%s %s %s filtered out spanning:%d minimum:%d\n' % (adj.contigs, 
	                                                                             adj.rearrangement,
	                                                                             adj_size,
	                                                                             support_total['spanning'],
	                                                                             min_support))
		
	if len(adj.homol_seq[0]) > max_homol_len:
	    if support_total['flanking'] < min_flanking:
		filtered_out = True
		if self.debug:
		    sys.stdout.write('%s %s %s homol_len:%d filtered out frags:%d minimum:%d\n' % (adj.contigs, 
		                                                                                   adj.rearrangement,
		                                                                                   adj_size,
		                                                                                   len(adj.homol_seq[0]),
		                                                                                   support_total['flanking'],
		                                                                                   min_flanking))	    
	
	if filtered_out and adj.novel_seq and adj.novel_seq != '-' and avg_tlen is not None:
	    spannable_fraction_tlen = 0.4
	    if len(adj.novel_seq) > spannable_fraction_tlen * avg_tlen:
		passed = False
		for tiling in support['tiling']:
		    if tiling:
			passed = True
			break
		filtered_out = True if not passed else False
		    
		if self.debug:
		    message = 'no_tiling_for_long_novel' if not passed else 'long_novel rescued by tiling'
		    sys.stdout.write('%s %s %s %s novel_seq:%s novel_len:%d %f tiling:%s spanning:%d frags:%d\n' % (message,
		                                                                           adj.contigs, 
		                                                                           adj.chroms, 
		                                                                           adj.breaks, 
		                                                                           adj.novel_seq,
		                                                                           len(adj.novel_seq), 
		                                                                           spannable_fraction_tlen * self.avg_tlen,
		                                                                           support['tiling'], 
		                                                                           support_total['spanning'],
		                                                                           support_total['flanking']))
	    	    	    	  
	if not normal:
	    adj.support_final = support_final
	else:
	    adj.support_final_normal = support_final

	return filtered_out
    
    def screen_by_coordinate(self, adjs, bad_bed_file):
	# creates bed file for all adjacencies
	adjs_bed_file = '%s/adjs.bed' % self.out_dir
	self.create_bed(adjs, adjs_bed_file)
	adjs_bed = BedTool(adjs_bed_file)
	
	# overlaps adjacencies with bad bed file, keeps breakpoints in Set
	bad_breaks = Set()
	regions = BedTool(bad_bed_file)
	overlaps = regions.intersect(adjs_bed)
	for olap in overlaps:
	    bad_breaks.add('%s:%s' % (olap[0], int(olap[1]) + 1))
	
	# screen out adjacencies
	bad_adj_indices = Set()
	for i in range(len(adjs)):
	    for j in (0,1):
		breakpt = '%s:%s' % (adjs[i].chroms[j], adjs[i].breaks[j])
		if breakpt in bad_breaks:
		    bad_adj_indices.add(i)
		    if self.debug:
			sys.stdout.write('%s %s:%s %s:%s (%s) overlaps repeat/segdup\n' % (adjs[i].contigs[0],
			                                                                      adjs[i].chroms[0],
			                                                                      adjs[i].breaks[0],
			                                                                      adjs[i].chroms[1],
			                                                                      adjs[i].breaks[1],
			                                                                      breakpt
			                                                                      )
			                 )
	for i in sorted(bad_adj_indices, reverse=True):
	    del adjs[i]
	    
    def create_bed(self, adjs, outfile):
	out = open(outfile, 'w')
	for adj in adjs:
	    out.write('%s\n' % adj.as_bed())
	out.close()
	  
    def filter_variants(self, min_support, min_support_normal=None, max_homol=None, tumor_bam=False, normal_bam=False):
	"""Filter out events that are believed to be false positive
	
	Filtering:
	1. Must not have 'N' in probe sequence
	2. Spanning reads must be >= minimum (default=4)
	3. If homology sequence > maximum (default=25), must have flanking read pair support
	
	Args:
	    min_spanning: (int) minimum number of spanning reads
	    min_flanking: (int) minimum number of flanking pairs when homlogous sequence > max_homol_len
	    max_homol_len: (int) maximum length of homology, over which flanking pairs have to be check
	"""
	for variant in self.variants:
	    for adj in variant.adjs:
		if 'N' in adj.probes[0].upper():
		    adj.filtered_out = True
		    if self.debug:
			sys.stdout.write('N_in_probe %s %s\n' % (adj.contigs, adj.probes[0]))
			
		if not adj.filtered_out and adj.novel_seq is not None and adj.novel_seq != 'NA' and adj.novel_seq != '-' and\
		     re.search('[^ATGC]', adj.novel_seq, re.IGNORECASE):
		    adj.filtered_out = True
		    if self.debug:
			sys.stdout.write('non AGTC in novel sequence %s %s\n' % (adj.contigs, adj.novel_seq))
			
		if not adj.filtered_out and (target_non_canonical(adj.chroms[0]) or target_non_canonical(adj.chroms[1])):
		    adj.filtered_out = True
		    if self.debug:
			sys.stdout.write('non canonical chromosome %s %s\n' % (adj.contigs, adj.chroms))
			
		if not adj.filtered_out and adj.homol_seq and max_homol is not None:
		    if len(adj.homol_seq[0]) > max_homol:
			adj.filtered_out = True
			if self.debug:
			    sys.stdout.write('homolgous sequence length %s too long (>%s)\n' % (len(adj.homol_seq[0]), max_homol))
			
		if not adj.filtered_out and (tumor_bam or normal_bam):
		    if tumor_bam:
			adj.filtered_out = self.filter_by_read_support(adj, min_support)
			
		    if normal_bam:
			adj.filtered_out_normal = self.filter_by_read_support(adj, min_support_normal, normal=True)
		    
			if adj.filtered_out is not None:
			    if not adj.filtered_out:
				adj.somatic = True
			    if adj.filtered_out_normal is not None and not adj.filtered_out_normal:
				adj.somatic = False
				
			elif adj.filtered_out_normal is not None:
			    if adj.filtered_out_normal:
				adj.somatic = True
							   
			if self.debug:
			    if not adj.somatic:
				sys.stdout.write('germline %s %s %s tumor:%s normal:%s\n' % (adj.contigs, 
				                                                             adj.chroms, 
				                                                             adj.breaks, 
				                                                             adj.support_total, 
				                                                             adj.support_total_normal))

	    for adj in variant.adjs:
		if adj.filtered_out:
		    variant.filtered_out = True
		    break
		
	    if normal_bam:
		for adj in variant.adjs:
		    if not adj.somatic:
			variant.somatic = False
			break
		    
def main(args, options):
    sv_finder = SVFinder(*args, 
                         genome=options.genome, 
                         index_dir=options.index_dir, 
                         num_procs=options.num_threads, 
                         is_transcriptome=options.is_transcriptome, 
                         annot_file=options.annot_file,
                         skip_simple_repeats=options.skip_simple_repeats,
                         cytobands_file=options.cytobands_file,
                         acen_buffer=options.acen_buffer,
                         debug=options.debug)
        
    adjs = sv_finder.find_adjs(min_ctg_cov=options.min_ctg_cov, 
                               max_size=options.max_size, 
                               min_size=options.min_size, 
                               ins_as_ins=options.ins_as_ins,
                               skip_acen=options.skip_acen,
                               check_alt_paths=options.check_alt_paths,
                               min_ctg_size=options.min_ctg_size,
                               bad_coords=options.bad_coords,
                               skip_contigs_file=options.skip_contigs)
    
    sv_finder.create_variants(adjs)
        
    if options.genome and options.index_dir and sv_finder.variants:
	sv_finder.screen_realigns(use_realigns=options.use_realigns)
	
    if options.r2c_bam_file or options.normal_bam:
	sv_finder.find_support(tumor_bam=options.r2c_bam_file, 
	                       min_overlap=options.min_overlap, 
	                       normal_bam=options.normal_bam, 
	                       min_overlap_normal=options.min_overlap_normal)
	    
    sv_finder.filter_variants(options.min_support, 
                              options.min_support_normal,
                              options.max_homol,
                              options.r2c_bam_file,
                              options.normal_bam)
    
	
    sv_finder.output(options.normal_bam is not None,
                     reference_url=options.reference_url,
                     assembly_url=options.assembly_url,
                     insertion_as_breakends=options.insertion_as_breakends,
                     )
    

if __name__ == '__main__':
    usage = "Usage: %prog c2g_bam aligner contig_fasta(indexed) genome_file(indexed) out_dir"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-T", "--is_transcriptome", dest="is_transcriptome", help="is sample transcriptome", action="store_true", default=False)
    parser.add_option("-a", "--annot_file", dest="annot_file", help="annotation file")
    parser.add_option("-b", "--r2c_bam", dest="r2c_bam_file", help="reads-to-contigs bam file")
    parser.add_option("-g", "--genome", dest="genome", help="genome")
    parser.add_option("-G", "--index_dir", dest="index_dir", help="genome index directory")
    parser.add_option("-t", "--num_threads", dest="num_threads", help="number of threads. Default:8", type='int', default=8)
    parser.add_option("--debug", dest="debug", help="debug mode", action="store_true", default=False)
    parser.add_option("--max_homol", dest="max_homol", help="maximum bases of microhomology. Default:25", type = "int", default=25)
    parser.add_option("--min_ctg_cov", dest="min_ctg_cov", help="minimum contig coverage. Default:0.95", type='float', default=0.95)
    parser.add_option("--min_support", dest="min_support", help="minimum read support. Default:2", type='int', default=2)
    parser.add_option("--min_support_normal", dest="min_support_normal", help="minimum read support for normal. Default:1", type='int', default=1)
    parser.add_option("--min_overlap", dest="min_overlap", help="minimum breakpoint overlap for identifying read support. Default:4", type='int', default=4)
    parser.add_option("--min_overlap_normal", dest="min_overlap_normal", help="minimum breakpoint overlap for identifying read support. Default:4", type='int', default=4)
    parser.add_option("--normal_bam", dest="normal_bam", help="reads-to-contigs bam file of match normal")
    parser.add_option("--reference_url", dest="reference_url", help="URL of reference")
    parser.add_option("--assembly_url", dest="assembly_url", help="URL of assembly file for VCF")
    parser.add_option("--skip_simple_repeats", dest="skip_simple_repeats", help="skip simple repeats", action="store_true", default=False)
    parser.add_option("--insertion_as_breakends", dest="insertion_as_breakends", help="outputs large insertion as breakends in VCF (Default will output as SV)", action="store_true", default=False)
    parser.add_option("--max_size", dest="max_size", help="maximum size of variant", type='int')
    parser.add_option("--min_size", dest="min_size", help="minimum size of variant", type='int')
    parser.add_option("--ins_as_ins", dest="ins_as_ins", help="keep small duplications as insertions", action="store_true", default=False)
    parser.add_option("--use_realigns", dest="use_realigns", help="use existing realignments", action="store_true", default=False)
    parser.add_option("--skip_acen", dest="skip_acen", help="skip acentromeric regions", action="store_true", default=False)
    parser.add_option("--cytobands", dest="cytobands_file", help="cytobands file")
    parser.add_option("--acen_buffer", dest="acen_buffer", help="buffer added to acen region. Default:100000", type="int", default=100000)
    parser.add_option("--check_alt_paths", dest="check_alt_paths", help="combine primary and secondary alignments to check for alternative path", action="store_true", default=False)
    parser.add_option("--min_ctg_size", dest="min_ctg_size", help="minimum contig size. Default:0 (no screening)", type="int", default=0)
    parser.add_option("--bad_coords", dest="bad_coords", help="BED file for coordinates to screen out e.g. segdups")
    parser.add_option("--skip_contigs", dest="skip_contigs", help="text file of contig names to skip")
    
    (options, args) = parser.parse_args()
    if len(args) == 5:
        main(args, options)


