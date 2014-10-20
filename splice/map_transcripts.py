from optparse import OptionParser
from itertools import groupby, chain
from pybedtools import create_interval_from_list, set_tempdir, BedTool
import pysam
import sys
import re
import os
from distutils.spawn import find_executable
from shared import gmap
from shared.annotate import overlap
from shared.alignment import reverse_complement
from sets import Set
from intspan import intspan
from SV.variant import Adjacency
from itd_finder import ITD_Finder
from fusion_finder import FusionFinder
from novel_splice_finder import NovelSpliceFinder
from shared.read_support import scan_all, fetch_support

class Transcript:
    def __init__(self, id, gene=None, strand=None):
	self.id = id
	self.gene = gene
	self.strand = strand
	self.exons = []
	
    def add_exon(self, exon):
	self.exons.append(exon)
	self.exons.sort(key=lambda e: int(e[0]))
	    
    def exon(self, num):
	assert type(num) is int, 'exon number %s not given in int' % num
	assert self.strand == '+' or self.strand == '-', 'transcript strand not valid: %s %s' % (self.id, self.strand)
	assert num >= 1 and num <= len(self.exons), 'exon number out of range:%s (1-%d)' % (num, len(self.exons))
	
	if self.strand == '+':
	    return self.exons[num - 1]
	else:
	    return self.exons[-1 * num]
	
    def num_exons(self):
	return len(self.exons)
    
    def length(self):
	total = 0
	for exon in self.exons:
	    total += (exon[1] - exon[0] + 1)
	return total
    
    def txStart(self):
	return self.exons[0][0]
    
    def txEnd(self):
	return self.exons[-1][1]
    
    def exon_num(self, index):
	"""Converts exon index to exon number
	Exon number is always starting from the transcription start site
	i.e. for positive transcripts, the first exon is exon 1
	     for negative transcripts, the last exon is exon 1
	Need this method because a lot of the splicing variant code just keep
	track of the index instead of actual exon number
	
	Args:
	    index: (int) index of exon in list
	Returns:
	    Exon number in int
	"""
	assert type(index) is int, 'exon index %s not given in int' % index
	assert self.strand == '+' or self.strand == '-', 'transcript strand not valid: %s %s' % (self.id, self.strand)
	assert index >= 0 and index < len(self.exons), 'exon index out of range:%s %d' % (index, len(self.exons))
	if self.strand == '+':
	    return index + 1
	else:
	    return len(self.exons) - index
    
class Event:
    # headers of tab-delimited output
    headers = ['event_type',
               'chrom1',
               'pos1',
               'orient1',
               'gene1',
               'transcript1',
               'exon1',
               'exon_bound1',
               'chrom2',
               'pos2',
               'orient2',
               'gene2',
               'transcript2',
               'exon2',
               'exon_bound2',
               'size',
               'novel_sequence',
               'homol_seq',
               'homol_coords',
               'homol_len',
               'contigs',
               'contig_breaks',
               'sense_fusion',
               "5'gene",
               "3'gene",
               'spanning_reads',
               ]
    
    @classmethod
    def output(cls, events, outdir):
	"""Output events

	Args:
	    events: (list) Adjacency
	    outdir: (str) absolute path of output directory
	Returns:
	    events will be output in file outdir/events.tsv
	"""
	out_file = '%s/events.tsv' % outdir
	out = open(out_file, 'w')
	out.write('%s\n' % '\t'.join(cls.headers))
	for event in events:
	    if event.rna_event:
		out_line = {
		    'fusion': cls.from_fusion,
		    'ITD': cls.from_single_locus,
		    'PTD': cls.from_fusion,
		    'ins': cls.from_single_locus,
		    'del': cls.from_single_locus,
		    'skipped_exon': cls.from_single_locus,
		    'novel_exon': cls.from_single_locus,
		    'novel_donor': cls.from_single_locus,
		    'novel_acceptor': cls.from_single_locus,
		    'novel_intron': cls.from_single_locus,
		    'retained_intron': cls.from_single_locus,
		    }[event.rna_event](event)
		if out_line:
		    out.write('%s\n' % out_line)
	out.close()
				
    @classmethod
    def from_fusion(cls, event):
	"""Generates output line for a fusion/PTD event
	
	Args:
	    event: (Adjacency) fusion or PTD event (coming from split alignment)
	Returns:
	    Tab-delimited line
	"""
	data = [event.rna_event]
	for values in zip(event.chroms, event.breaks, event.orients, event.genes, event.transcripts, event.exons, event.exon_bound):
	    data.extend(values)
	    
	# size not applicable to fusion
	data.append('-')
	if hasattr(event, 'novel_seq') and event.novel_seq is not None:
	    data.append(event.novel_seq)
	else:
	    data.append('-')
	    
	# homol_seq and coords
	homol_seq = None
	if type(event.homol_seq) is str:
	    homol_seq = event.homol_seq
	elif type(event.homol_seq) is list or type(event.homol_seq) is tuple:
	    homol_seq = event.homol_seq[0]
	if homol_seq:
	    data.append(homol_seq)
	
	    if homol_seq != '-':
		data.append('-'.join(map(str, event.homol_coords)))
		data.append(event.homol_coords[1] - event.homol_coords[0] + 1)
	    else:
		data.append('-')
		data.append('-')
	else:
	    data.append('-')
	    data.append('-')
	    data.append('-')
	
	data.append(','.join(event.contigs))
	data.append(cls.to_string(event.contig_breaks))
	
	data.append(event.is_sense)
	data.append(event.gene5)
	data.append(event.gene3)
	
	if not event.support['spanning']:
	    data.append('-')
	else:
	    data.append(event.support['spanning'][0])
	
	return '\t'.join(map(str, data))
    
    @classmethod
    def from_single_locus(cls, event):
	"""Generates output line for an event from a single alignment
	
	Args:
	    event: (Adjacency) indel, ITD, splicing event (coming from single alignment)
	Returns:
	    Tab-delimited line
	"""
	chroms = (event.chroms[0], event.chroms[0])
	orients = ('L', 'R')
	genes = (event.genes[0], event.genes[0])
	transcripts = (event.transcripts[0], event.transcripts[0])
	if event.exons:
	    if len(event.exons) == 2:
		exons = event.exons
	    else:
		exons = (event.exons[0], event.exons[0])
	# novel exons
	else:
	    exons = ('-', '-')
	    
	# exon_bound
	exon_bound = ('-', '-')
	
	data = [event.rna_event]
	for values in zip(chroms, event.breaks, orients, genes, transcripts, exons, exon_bound):
	    data.extend(values)
	
	if hasattr(event, 'size') and event.size is not None:
	    data.append(event.size)
	else:
	    data.append('-')
	if hasattr(event, 'novel_seq') and event.novel_seq is not None:
	    data.append(event.novel_seq)
	else:
	    data.append('-')
	    
	# homol_seq and coords
	data.append('-')
	data.append('-')
	data.append('-')
		
	data.append(','.join(event.contigs))
	data.append(cls.to_string(event.contig_breaks))
	
	# sense fusion
	data.append('-')
	data.append('-')
	data.append('-')
	
	if not event.support['spanning']:
	    print 'failed_spanning', event.rna_event, event.contigs, event.support
	    data.append('-')
	else:
	    data.append(event.support['spanning'][0])

	return '\t'.join(map(str, data))
	
    @classmethod
    def to_string(cls, value):
	"""Convert value of data types other than string usually used
	in Adjacency attributes to string for print
	
	Args:
	    value: can be
	           1. None
		   2. simple list/tuple
		   3. list/tuple of list/tuple
	Returns:
	    string representation of value
	    ';' used to separate items in top-level list/tuple
	    ',' used to separate items in second-level list/tuple
	"""
	if value is None:
	    return '-'
	elif type(value) is list or type(value) is tuple:
	    items = []
	    for item in value:
		if item is None:
		    items.append('-')
		elif type(item) is list or type(item) is tuple:
		    if item:
			items.append(','.join(map(str, item)))
		else:
		    items.append(str(item))
		    
	    if items:
		return ';'.join(items)
	    else:
		return '-'
	else:
	    return str(value)
	
    @classmethod
    def screen(cls, events, outdir, aligner, align_info=None, debug=False):
	"""Screen events identified and filter out bad ones
	
	Right now it just screens out fusion whose probe sequence can align to one single location
	
	Args:
	    events: (list) Events
	    outdir: (str) absolute path of output directory, for storing re-alignment results
	    aligner: (str) aligner name (gmap)
	    align_info: (dict) 'genome', 'index_dir', 'num_procs'
	    debug: (boolean) output debug info e.g. reason for screening out event
	"""
	fusions = [e for e in events if e.rna_event == 'fusion']
	bad_contigs = FusionFinder.screen_realigns(fusions, outdir, aligner, align_info, debug=debug)
	
	bad_event_indices = []
	# remove any event that involve contigs that failed screening as mapping is not reliable
	for e in reversed(range(len(events))):
	    for contig in events[e].contigs:
		if contig in bad_contigs and not e in bad_event_indices:
		    bad_event_indices.append(e)
		    break
		
	for e in bad_event_indices:
	    del events[e]
	
class Mapping:
    """Mapping per alignment"""
    def __init__(self, contig, align_blocks, transcripts=[]):
	self.contig = contig
	self.transcripts = transcripts
	self.genes = list(Set([txt.gene for txt in self.transcripts]))
	self.align_blocks = align_blocks
	# coverage for each transcript in self.transcripts
	self.coverages = []
	
    def overlap(self):
	"""Overlaps alignment spans with exon-spans of each matching transcripts
	Upates self.coverages
	"""
	align_span = self.create_span(self.align_blocks)
	
	for transcript in self.transcripts:
	    exon_span = self.create_span(transcript.exons)
	    olap = exon_span.intersection(align_span)
	    self.coverages.append(float(len(olap)) / float(len(exon_span)))
	    
    @classmethod
    def create_span(cls, blocks):
	"""Creates intspan for each block
	Used by self.overlap()
	"""
	span = None
	for block in blocks:
	    try:
		span = span.union(intspan('%s-%s' % (block[0], block[1])))
	    except:
		span = intspan('%s-%s' % (block[0], block[1]))
		
	return span
    
    @classmethod
    def header(cls):
	return '\t'.join(['contig',
	                  'gene',
	                  'transcript',
	                  'coverage'
	                  ])
	    
    def as_tab(self):
	"""Generates tab-delimited line for each mapping"""
	data = []
	data.append(self.contig)
	if self.transcripts:
	    data.append(','.join([gene for gene in Set([txt.gene for txt in self.transcripts])]))
	    data.append(','.join([txt for txt in Set([txt.id for txt in self.transcripts])]))
	else:
	    data.append('-')
	    data.append('-')
	    
	if self.coverages:
	    data.append(','.join(['%.2f' % olap for olap in self.coverages]))
	else:
	    data.append('-')
	
	return '\t'.join(map(str, data))
	
    @classmethod
    def pick_best(self, mappings, align, debug=False):
	"""Selects best mapping among transcripts"""
	scores = {}
	metrics = {}
	for transcript, matches in mappings:
	    metric = {'score': 0,
	              'from_edge': 0,
	              'txt_size': 0
	              }
	    # points are scored for matching exon boundaries
	    score = 0
	    for i in range(len(matches)):
		if matches[i] is None:
		    continue
		
		if i == 0:			    
		    if matches[i][0][1][0] == '=':
			score += 5
		    elif matches[i][0][1][0] == '>':
			score += 2
			
		elif i == len(matches) - 1:
		    if matches[i][-1][1][1] == '=':
			score += 5
		    elif matches[i][-1][1][1] == '<':
			score += 2
			
		if matches[i][0][1][1] == '=':
		    score += 4
		if matches[i][-1][1][1] == '=':
		    score += 4	
			
	    # points are deducted for events
	    penalty = 0
	    for i in range(len(matches)):
		# block doesn't match any exon
		if matches[i] is None:
		    penalty += 2
		    continue
		    
		# if one block is mapped to >1 exon
		if len(matches[i]) > 1:
		    penalty += 1
		
		# if consecutive exons are not mapped to consecutive blocks
		if i < len(matches) - 1 and matches[i + 1] is not None:
		    if matches[i + 1][0][0] != matches[i][-1][0] + 1:
			#print 'penalty', matches[i][-1][0], matches[i + 1][0][0]
			penalty += 1
			
	    metric['score'] = score - penalty
	    
	    # if the first or last block doesn't have matches, won't be able to calculate distance from edges
	    # in that case, will set the distance to 'very big'
	    if matches[0] is None or matches[-1] is None:
		metric['from_edge'] = 10000
	    else:
		if transcript.strand == '+':
		    start_exon_num = matches[0][0][0] + 1
		    end_exon_num = matches[-1][-1][0] + 1
		else:
		    start_exon_num = transcript.num_exons() - matches[0][0][0]
		    end_exon_num = transcript.num_exons() - matches[-1][-1][0]
		    
		metric['from_edge'] = align.tstart - transcript.exon(start_exon_num)[0] + align.tend - transcript.exon(end_exon_num)[1]
		    
	    metric['txt_size'] = transcript.length()
	    metrics[transcript] = metric
	    
	    if debug:
		sys.stdout.write("mapping %s %s %s %s %s %s\n" % (align.query, transcript.id, transcript.gene, score, penalty, metric))
	    
	transcripts_sorted = sorted(metrics.keys(), key = lambda txt: (-1 * metrics[txt]['score'], metrics[txt]['from_edge'], metrics[txt]['txt_size']))
	for t in transcripts_sorted:
	    print 'sorted', t.id, metrics[t]
	    	    
	best_transcript = transcripts_sorted[0]
	best_matches = [mapping[1] for mapping in mappings if mapping[0] == best_transcript]
	best_mapping = Mapping(align.query,
	                       align.blocks,
                               [transcripts_sorted[0]],
                               )
	best_mapping.overlap()
	
	return best_mapping
	    	
    @classmethod
    def group(cls, all_mappings):
	"""Group mappings by gene"""
	gene_mappings = []
	for gene, group in groupby(all_mappings, lambda m: m.genes[0]):
	    mappings = list(group)
	    contigs = ','.join([mapping.contig for mapping in mappings])
	    transcripts = [mapping.transcripts for mapping in mappings]
	    align_blocks = [mapping.align_blocks for mapping in mappings]
	    
	    align_blocks = None
	    
	    for mapping in mappings:
		try:
		    align_blocks = align_blocks.union(cls.create_span(mapping.align_blocks))
		except:
		    align_blocks = cls.create_span(mapping.align_blocks)		
	    
	    gene_mappings.append(Mapping(contigs,
	                                 align_blocks.ranges(),
	                                 list(Set(chain(*transcripts))),
	                                 )
	                         )
	    	    
	[mapping.overlap() for mapping in gene_mappings]
	return gene_mappings
    
    @classmethod
    def output(cls, mappings, outfile):
	"""Output mappings into output file"""
	out = open(outfile, 'w')
	out.write('%s\n' % cls.header())
	for mapping in mappings:
	    out.write('%s\n' % mapping.as_tab())
	out.close()

class ExonMapper:
    def __init__(self, bam_file, aligner, contigs_fasta_file, annotation_file, ref_fasta_file, outdir, 
                 itd_min_len=None, itd_min_pid=None, itd_max_apart=None, debug=False):
        self.bam = pysam.Samfile(bam_file, 'rb')
	self.contigs_fasta_file = contigs_fasta_file
	self.contigs_fasta = pysam.Fastafile(contigs_fasta_file)
        self.ref_fasta = pysam.Fastafile(ref_fasta_file)
	self.annot = pysam.Tabixfile(annotation_file, parser=pysam.asGTF())
        self.aligner = aligner
        self.outdir = outdir
	self.debug = debug
        
        self.blocks_bed = '%s/blocks.bed' % outdir
        self.overlaps_bed = '%s/blocks_olap.bed' % outdir
        self.aligns = {}        

        self.annotations_file = annotation_file
	
	self.mappings = []
	self.events = []
	
	# initialize ITD conditions
	self.itd_conditions = {'min_len': itd_min_len,
	                       'min_pid': itd_min_pid,
	                       'max_apart': itd_max_apart
	                       }
					
    def map_contigs_to_transcripts(self):
	"""Maps contig alignments to transcripts, discovering variants at the same time"""
	# extract all transcripts info in dictionary
	transcripts = self.extract_transcripts()
	
	aligns = []
	for contig, group in groupby(self.bam.fetch(until_eof=True), lambda x: x.qname):
	    sys.stdout.write('analyzing %s\n' % contig)
            alns = list(group)	    
	    aligns = self.extract_aligns(alns)
	    if aligns is None:
		sys.stdout.write('no valid alignment: %s\n' % contig)
		continue
	    	    
	    # for finding microhomolgy sequence and generating probe in fusion
	    contig_seq = self.contigs_fasta.fetch(contig)
	    
	    chimera = True if len(aligns) > 1 else False
	    chimera_block_matches = []
	    for align in aligns:
		if align is None:
		    sys.stdout.write('bad alignment: %s\n' % contig)
		    continue
		
		if re.search('[._Mm]', align.target):
		    sys.stdout.write('skip target:%s %s\n' % (contig, align.target))
		    continue
		
		if self.debug:
		    sys.stdout.write('contig:%s genome_blocks:%s contig_blocks:%s\n' % (align.query, 
		                                                                        align.blocks, 
		                                                                        align.query_blocks
		                                                                        ))
		
		# entire contig align within single exon or intron
		within_intron = []
		within_exon = []
		
		transcripts_mapped = Set()
		events = []
		# each gtf record corresponds to a feature
		for gtf in self.annot.fetch(align.target, align.tstart, align.tend):	
		    # contigs within single intron or exon
		    if not chimera and align.tstart >= gtf.start and align.tend <= gtf.end:
			match = self.match_exon((align.tstart, align.tend), (gtf.start, gtf.end)) 
			
			if gtf.feature == 'intron':	
			    within_intron.append((gtf, match))
			    
			elif gtf.feature == 'exon':
			    within_exon.append((gtf, match))
			  
		    # collect all the transcripts that have exon overlapping alignment
		    else:
			if gtf.feature == 'exon':
			    transcripts_mapped.add(gtf.transcript_id)
			    		
		if transcripts_mapped:
		    mappings = []
		    # key = transcript name, value = "matches"
		    all_block_matches = {}
		    # Transcript objects that are fully matched
		    full_matched_transcripts = []
		    for txt in transcripts_mapped:
			block_matches = self.map_exons(align.blocks, transcripts[txt].exons)
			all_block_matches[txt] = block_matches
			mappings.append((transcripts[txt], block_matches))
			
			if not chimera and self.is_full_match(block_matches):
			    full_matched_transcripts.append(transcripts[txt])
			    
		    if not full_matched_transcripts:	
			best_mapping = Mapping.pick_best(mappings, align, debug=True)
			self.mappings.append(best_mapping)
			
			# find events only for best transcript
			best_transcript = best_mapping.transcripts[0]
			events = self.find_events({best_transcript.id:all_block_matches[best_transcript.id]}, 
			                          align, 
			                          {best_transcript.id:best_transcript})
			#events = self.find_events(all_block_matches, align, transcripts)
			if events:
			    self.events.extend(events)
			elif self.debug:
			    sys.stdout.write('%s - partial but no events\n' % align.query)		    
			    			    
		    if chimera:
			chimera_block_matches.append(all_block_matches)
		    
		elif not chimera:
		    if within_exon:
			sys.stdout.write("contig mapped within single exon: %s %s:%s-%s %s\n" % (contig, 
			                                                                         align.target, 
			                                                                         align.tstart, 
			                                                                         align.tend, 
			                                                                         within_exon[0]
			                                                                         ))
		    
		    elif within_intron:
			sys.stdout.write("contig mapped within single intron: %s %s:%s-%s %s\n" % (contig, 
			                                                                           align.target, 
			                                                                           align.tstart, 
			                                                                           align.tend, 
			                                                                           within_intron[0]
			                                                                           ))
		
	    # split aligns, try to find gene fusion
	    if chimera and chimera_block_matches:
		if len(chimera_block_matches) == len(aligns):
		    fusion = FusionFinder.find_chimera(chimera_block_matches, transcripts, aligns, contig_seq)
		    if fusion:
			homol_seq, homol_coords = None, None
			if self.aligner.lower() == 'gmap':
			    homol_seq, homol_coords = gmap.find_microhomology(alns[0], contig_seq)
			if homol_seq is not None:
			    fusion.homol_seq = homol_seq
			    fusion.homol_coords = homol_coords
			self.events.append(fusion)
		
	    ## pick best matches
	    #if len(mappings) > 1:
		#best_mapping = Mapping.pick_best(mappings, align, debug=True)
		#self.mappings.append(best_mapping)
					
    def map_exons(self, blocks, exons):
	"""Maps alignment blocks to exons
	
	Crucial in mapping contigs to transcripts
	
	Args:
	    blocks: (List) of 2-member list (start and end of alignment block)
	    exons: (List) of 2-member list (start and end of exon)
	Returns:
	    List of list of tuples, where
	         each item of the top-level list corresponds to each contig block
		 each item of the second-level list corresponds to each exon match
		 each exon-match tuple contains exon_index, 2-character matching string
		 e.g. [ [ (0, '>='), (1, '<=') ], [ (2, '==') ], [ (3, '=<') ], [ (3, '>=') ] ]
		 this says that the alignment has 4 blocks,
		 first block matches to exon 0 and 1, find a retained-intron,
		 second block matches perfectly to exon 2
		 third and fourth blocks matching to exon 3, possible a deletion or novel_intron
	"""
	result = []
	for b in range(len(blocks)):
	    block_matches = []
	    
	    for e in range(len(exons)):
		block_match = self.match_exon(blocks[b], exons[e])
		if block_match != '':
		    block_matches.append((e, block_match))
		    
	    if not block_matches:
		block_matches = None
	    result.append(block_matches)
		
	return result
		    		    
    def extract_transcripts(self):
	"""Extracts all exon info into transcript objects
	
	Requires annotation file passed to object
	Uses PyBedTool for parsing
	
	Returns:
	    List of Transcripts with exon info, strand
	"""
	transcripts = {}
	for feature in BedTool(self.annotations_file):
	    if feature[2] == 'exon':
		exon = (int(feature.start) + 1, int(feature.stop))
		exon_num = int(feature.attrs['exon_number'])
		transcript_id = feature.attrs['transcript_id']
		gene = feature.attrs['gene_name']
		strand = feature.strand
				
		try:
		    transcript = transcripts[transcript_id]
		except:
		    transcript = Transcript(transcript_id, gene=gene, strand=strand)
		    transcripts[transcript_id] = transcript
		    
		transcript.add_exon(exon)
		
	return transcripts
        
    def extract_novel_seq(self, adj):
	"""Extracts novel sequence in adjacency, should be a method of Adjacency"""
	contig_breaks = adj.contig_breaks[0]
	start, end = (contig_breaks[0], contig_breaks[1]) if contig_breaks[0] < contig_breaks[1] else (contig_breaks[1], contig_breaks[0])
	novel_seq = self.contigs_fasta.fetch(adj.contigs[0], start, end - 1)
	return novel_seq
	    
    def find_events(self, matches_by_transcript, align, transcripts):
	"""Find events from single alignment
	
	Wrapper for finding events within a single alignment
	Maybe a read-through fusion, calls FusionFinder.find_read_through()
	Maybe splicing or indels, call NovelSpliceFinder.find_novel_junctions()
	Will take results in dictionary and convert them to Adjacencies
	
	Args:
	    matches_by_transcript: (list) dictionaries where 
	                                      key=transcript name, 
	                                      value=[match1, match2, ...] where
						    match1 = matches of each alignment block
							     i.e.
							     [(exon_id, '=='), (exon_id, '==')] 
							     or None if no_match
	    align: (Alignment) alignment object
	    transcripts: (dictionary) key=transcript_id value: Transcript
	Returns:
	    List of Adjacencies
	"""
	events = []
	
	# for detecting whether there is a read-through fusion
	genes = Set([transcripts[txt].gene for txt in matches_by_transcript.keys()])		
	if len(genes) > 1:
	    fusion = FusionFinder.find_read_through(matches_by_transcript, transcripts, align)
	    if fusion is not None:
		events.append(fusion)
	    	    	
	# events within a gene
	local_events = NovelSpliceFinder.find_novel_junctions(matches_by_transcript, align, transcripts, self.ref_fasta)
	if local_events:
	    for event in local_events:
		adj = Adjacency((align.target,), event['pos'], '-', contig=align.query)
		if event['event'] in ('ins', 'del', 'dup', 'inv'):
		    adj.rearrangement = event['event']
		adj.rna_event = event['event']
		adj.genes = (list(genes)[0],)
		adj.transcripts = (event['transcript'][0],)
		adj.size = event['size']
		
		# converts exon index to exon number (which takes transcript strand into account)
		exons = event['exons'][0]
		#if type(exons) is not list and type(exons) is not tuple:
		    #print 'check', align.query, event, exons
		    #continue
		adj.exons = map(transcripts[event['transcript'][0]].exon_num, exons)
		
		adj.contig_breaks = event['contig_breaks']
		
		if adj.rearrangement == 'ins' or adj.rearrangement == 'dup':
		    novel_seq = self.extract_novel_seq(adj)
		    adj.novel_seq = novel_seq if align.strand == '+' else reverse_complement(novel_seq)
		    
		    if len(novel_seq) >= self.itd_conditions['min_len']:
			ITD_Finder.detect_itd(adj, align, self.contigs_fasta.fetch(adj.contigs[0]), self.outdir, 
			                      self.itd_conditions['min_len'],
			                      self.itd_conditions['max_apart'],
			                      self.itd_conditions['min_pid'],
			                      debug=self.debug
			                      )
		    			
		    adj.size = len(adj.novel_seq)
		    
		elif adj.rearrangement == 'del':
		    adj.size = adj.breaks[1] - adj.breaks[0] - 1
		
		events.append(adj)
			    
	return events
		            
    def extract_aligns(self, alns):
	"""Generates Alignments objects of chimeric and single alignments

	Args:
	    alns: (list) Pysam.AlignedRead
	Returns:
	    List of Alignments that are either chimeras or single alignments
	"""
	try:
	    chimeric_aligns = {
		'gmap': gmap.find_chimera,
		}[self.aligner](alns, self.bam)
	    if chimeric_aligns:
		return chimeric_aligns
	    else:
		return [{
		    'gmap': gmap.find_single_unique,
		    }[self.aligner](alns, self.bam)]
	except:
	    sys.exit("can't convert \"%s\" alignments - abort" % self.aligner)
        
    def match_exon(self, block, exon):
	"""Match an alignment block to an exon
	
	Args:
	    block: (Tuple/List) start and end of an alignment block
	    exon: (Tuple/List) start and end of an exon
	Returns:
	    2 character string, each character the result of each coordinate
	    '=': block coordinate = exon coordinate
	    '<': block coordinate < exon coordinate
	    '>': block coordinate > exon coordinate
	"""
        assert len(block) == 2 and len(block) == len(exon), 'unmatched number of block(%d) and exon(%d)' % (len(block), len(exon))
        assert type(block[0]) is int and type(block[1]) is int and type(exon[0]) is int and type(exon[1]) is int,\
        'type of block and exon must be int'
        result = ''
        
	if min(block[1], exon[1]) - max(block[0], exon[0]) > 0:
	    for i in range(0, 2):
		if block[i] == exon[i]:
		    result += '='
		elif block[i] > exon[i]:
		    result += '>'
		else:
		    result += '<'
        
        return result
    
    def is_full_match(self, block_matches):
	"""Determines if a contig fully covers a transcript

	A 'full' match is when every exon boundary matches, except the terminal boundaries
	
	Args:
	    block_matches: List of list of tuples, where
			   each item of the top-level list corresponds to each contig block
			   each item of the second-level list corresponds to each exon match
			   each exon-match tuple contains exon_index, 2-character matching string
			   e.g. [ [ (0, '>='), (1, '<=') ], [ (2, '==') ], [ (3, '=<') ], [ (3, '>=') ] ]
			   this says that the alignment has 4 blocks,
			   first block matches to exon 0 and 1, find a retained-intron,
			   second block matches perfectly to exon 2
			   third and fourth blocks matching to exon 3, possible a deletion or novel_intron

	Returns:
	    True if full, False if partial
	"""
	if None in block_matches:
	    return False
	
	if len(block_matches) == 1:
	    if len(block_matches[0]) == 1:
		if block_matches[0][0][1] == '==' or block_matches[0][0][1] == '>=' or block_matches[0][0][1] == '=<':
		    return True
	    
	    return False
	
	# if a block is mapped to >1 exon
	if [m for m in block_matches if len(m) > 1]:
	    return False
	
	exons = [m[0][0] for m in block_matches]
	
	if block_matches[0][0][1][1] == '=' and\
	   block_matches[-1][0][1][0] == '=' and\
	   len([m for m in block_matches[1:-1] if m[0][1] == '==']) == len(block_matches) - 2 and\
	   (len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(block_matches) - 1 or\
	    len([(a, b) for a, b in zip(exons, exons[1:]) if b == a - 1]) == len(block_matches) - 1):
	    return True
		
	return False
    
    def find_support(self, r2c_bam_file=None, min_overlap=None, multi_mapped=False, perfect=True, num_procs=1):
	"""Extracts read support from reads-to-contigs BAM
	
	Assumes reads-to-contigs is NOT multi-mapped 
	(so we add the support of different contigs of the same event together)
	
	Args:
	    bam_file: (str) Path of reads-to-contigs BAM
	"""
	coords = {}
	for event in self.events:
	    for i in range(len(event.contigs)):
		contig = event.contigs[i]
		print 'check', event.rna_event, contig, event.contigs, i, event.contig_breaks
		span = event.contig_breaks[i][0], event.contig_breaks[i][1]
		try:
		    coords[contig].append(span)
		except:
		    coords[contig] = [span]
	
	if coords:	    
	    avg_tlen = None
	    tlens = []
	    # if fewer than 10000 adjs, use 'fetch' of Pysam
	    if len(coords) < 1000 or num_procs == 1:
		support, tlens = fetch_support(coords, r2c_bam_file, self.contigs_fasta, overlap_buffer=min_overlap, perfect=perfect, debug=self.debug)
	    # otherwise use multi-process parsing of bam file
	    else:
		support, tlens = scan_all(coords, r2c_bam_file, self.contigs_fasta_file, num_procs, overlap_buffer=min_overlap, perfect=perfect, debug=self.debug)
		    
	    print 'total tlens', len(tlens)
		
	    for event in self.events:
		for i in range(len(event.contigs)):
		    contig = event.contigs[i]
		    span = event.contig_breaks[i][0], event.contig_breaks[i][1]
		    coord = '%s-%s' % (span[0], span[1])
					
		    if support.has_key(contig) and support[contig].has_key(coord):
			event.support['spanning'].append(support[contig][coord][0])
			
		if not multi_mapped:
		    event.sum_support()	
		    
		print event.support
		     
	    if tlens:
		avg_tlen = float(sum(tlens)) / len(tlens)
		print 'avg tlen', avg_tlen, sample_type
		    	
	if self.debug:
	    for event in self.events:
		print 'support', event.support, event.support_total, event.support_final


def main(args, options):
    outdir = args[-1]
    # check executables
    required_binaries = ('gmap', 'blastn')
    for binary in required_binaries:
	which = find_executable(binary)
	if not which:
	    sys.exit('"%s" not in PATH - abort' % binary)
    
    # find events
    em = ExonMapper(*args, 
                    itd_min_len=options.itd_min_len,
                    itd_min_pid=options.itd_min_pid,
                    itd_max_apart=options.itd_max_apart,
                    debug=options.debug)
    em.map_contigs_to_transcripts()
    
    Event.screen(em.events, outdir, em.aligner, {'genome': options.genome,
                                                 'index_dir': options.index_dir,
                                                 'num_procs': options.num_threads,
                                                 },
                 debug=options.debug)
    
    if options.r2c_bam_file:
	em.find_support(options.r2c_bam_file, options.min_overlap, options.multimapped, num_procs=options.num_threads)

    Mapping.output(em.mappings, '%s/contig_mappings.tsv' % outdir)
    gene_mappings = Mapping.group(em.mappings)
    Mapping.output(gene_mappings, '%s/gene_mappings.tsv' % outdir)
    
    Event.output(em.events, outdir)
    
if __name__ == '__main__':
    usage = "Usage: %prog c2g_bam aligner contigs_fasta annotation_file genome_file(indexed) out_dir"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-b", "--r2c_bam", dest="r2c_bam_file", help="reads-to-contigs bam file")
    parser.add_option("-t", "--num_threads", dest="num_threads", help="number of threads. Default:8", type='int', default=8) 
    parser.add_option("-g", "--genome", dest="genome", help="genome")
    parser.add_option("-G", "--index_dir", dest="index_dir", help="genome index directory")
    parser.add_option("--junctions", dest="junctions", help="output junctions", action="store_true", default=False)
    parser.add_option("--itd_min_len", dest="itd_min_len", help="minimum ITD length. Default: 10", default=10, type=int)
    parser.add_option("--itd_min_pid", dest="itd_min_pid", help="minimum ITD percentage of identity. Default: 0.95", default=0.95, type=float)
    parser.add_option("--itd_max_apart", dest="itd_max_apart", help="maximum distance apart of ITD. Default: 10", default=10, type=int)
    parser.add_option("--multimapped", dest="multimapped", help="reads-to-contigs alignment is multi-mapped", action="store_true", default=False)
    parser.add_option("--min_overlap", dest="min_overlap", help="minimum breakpoint overlap for identifying read support. Default:4", type='int', default=4)
    parser.add_option("--debug", dest="debug", help="debug mode", action="store_true", default=False)
    
    (options, args) = parser.parse_args()
    if len(args) == 6:
        main(args, options)     