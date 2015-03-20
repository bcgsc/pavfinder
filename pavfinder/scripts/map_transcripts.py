from optparse import OptionParser
from itertools import groupby
import pysam
import sys
import re
import os
import gzip
from sets import Set
from distutils.spawn import find_executable
from pavfinder import __version__
from pavfinder.shared import gmap
from pavfinder.shared.alignment import Alignment, reverse_complement
from pavfinder.shared.adjacency import Adjacency
from pavfinder.splice.event import Event
from pavfinder.splice.mapping import Mapping
from pavfinder.splice.transcript import Transcript
from pavfinder.splice import itd_finder
from pavfinder.splice import fusion_finder
from pavfinder.splice import novel_splice_finder
from pavfinder.shared.read_support import scan_all, fetch_support, expand_contig_breaks

class ExonMapper:
    def __init__(self, bam_file, aligner, contigs_fasta_file, annotation_file, ref_fasta_file, outdir, 
                 itd_min_len=None, itd_min_pid=None, itd_max_apart=None, 
                 exon_bound_fusion_only=False, coding_fusion_only=False, sense_fusion_only=False,
                 debug=False):
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

        self.annotation_file = annotation_file
	
	self.mappings = []
	self.events = []
	
	# initialize ITD conditions
	self.itd_conditions = {'min_len': itd_min_len,
	                       'min_pid': itd_min_pid,
	                       'max_apart': itd_max_apart
	                       }
	
	self.fusion_conditions = {'exon_bound_only': exon_bound_fusion_only,
	                          'coding_only': coding_fusion_only,
	                          'sense_only': sense_fusion_only}

    def map_contigs_to_transcripts(self):
	"""Maps contig alignments to transcripts, discovering variants at the same time"""
	# extract all transcripts info in dictionary
	transcripts = Transcript.extract_transcripts(self.annotation_file)

	chimeras = {}
	for aln in self.bam.fetch(until_eof=True):
	    if aln.is_unmapped:
		continue

	    contig = aln.query_name
	    contig_seq = self.contigs_fasta.fetch(contig)

	    align = Alignment.from_alignedRead(aln, self.bam)
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
		# collect all the transcripts that have exon overlapping alignment
		if gtf.feature == 'exon' or gmap.is_chimera(aln):
		    transcripts_mapped.add(gtf.transcript_id)
		# contigs within single intron
		elif gtf.feature == 'intron' and\
	             align.tstart >= gtf.start and align.tend <= gtf.end:
		    match = self.match_exon((align.tstart, align.tend), (gtf.start, gtf.end))
		    within_intron.append((gtf, match))

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

		    if self.is_full_match(block_matches):
			full_matched_transcripts.append(transcripts[txt])

		# report mapping
		best_mapping = Mapping.pick_best(mappings, align, debug=self.debug)
		self.mappings.append(best_mapping)

		if not full_matched_transcripts:
		    # find events only for best transcript
		    best_transcript = best_mapping.transcripts[0]
		    events = self.find_events({best_transcript.id:all_block_matches[best_transcript.id]},
		                              align,
		                              {best_transcript.id:best_transcript})
		    for event in events:
			event.contig_sizes.append(len(contig_seq))
		    if events:
			self.events.extend(events)
		    elif self.debug:
			sys.stdout.write('%s - partial but no events\n' % align.query)

		if gmap.is_chimera(aln):
		    try:
			chimeras[contig].append((align, all_block_matches))
		    except:
			chimeras[contig] = [(align, all_block_matches)]

	    else:
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
	if chimeras:
	    for contig in sorted(chimeras.keys()):
		# in case a chimeric alignment is entirely within intron
		# and introns are not features in given gtf, it will be not captured
		# and there will not be two chimeras
		if len(chimeras[contig]) < 2:
		    continue

		contig_seq = self.contigs_fasta.fetch(contig)
		aligns = [chimeras[contig][0][0], chimeras[contig][1][0]]
		chimera_block_matches = [chimeras[contig][0][1], chimeras[contig][1][1]]
		# order is required to be correct
		if aligns[0].qstart > aligns[1].qstart:
		    aligns.reverse()
		    chimera_block_matches.reverse()

		fusion = fusion_finder.find_chimera(chimera_block_matches, transcripts, aligns, contig_seq,
		                                   exon_bound_only=self.fusion_conditions['exon_bound_only'],
		                                   coding_only=self.fusion_conditions['coding_only'],
		                                   sense_only=self.fusion_conditions['sense_only'])
		if fusion:
		    homol_seq, homol_coords = None, None
		    if self.aligner.lower() == 'gmap':
			homol_seq, homol_coords = gmap.find_microhomology(aligns[0].sam, contig_seq)
		    if homol_seq is not None:
			fusion.homol_seq.append(homol_seq)
			fusion.homol_coords.append(homol_coords)
		    fusion.contig_sizes.append(len(contig_seq))
		    self.events.append(fusion)

	# expand contig span
	for event in self.events:
	    # if contig_support_span is not defined (it can be pre-defined in ITD)
	    # then set it
	    if not event.contig_support_span:
		event.contig_support_span = event.contig_breaks
		if event.rearrangement == 'ins' or event.rearrangement == 'dup':
		    expanded_contig_breaks = expand_contig_breaks(event.chroms[0],
			                                          event.breaks,
			                                          event.contigs[0],
			                                          [event.contig_breaks[0][0] + 1, event.contig_breaks[0][1] - 1],
			                                          event.rearrangement,
			                                          self.ref_fasta,
			                                          self.contigs_fasta,
			                                          self.debug)
		    if expanded_contig_breaks is not None:
			event.contig_support_span = [(expanded_contig_breaks[0] - 1, expanded_contig_breaks[1] + 1)]

    def map_contigs_to_transcripts_old(self):
	"""Maps contig alignments to transcripts, discovering variants at the same time"""
	# extract all transcripts info in dictionary
	transcripts = Transcript.extract_transcripts(self.annotation_file)
	
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
		    # collect all the transcripts that have exon overlapping alignment
		    if gtf.feature == 'exon':
			transcripts_mapped.add(gtf.transcript_id)
		    # contigs within single intron 
		    elif gtf.feature == 'intron' and\
		         not chimera and\
		         align.tstart >= gtf.start and align.tend <= gtf.end:
			match = self.match_exon((align.tstart, align.tend), (gtf.start, gtf.end)) 
			within_intron.append((gtf, match))		
			    	
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
			    
		    # report mapping
		    best_mapping = Mapping.pick_best(mappings, align, debug=self.debug)
		    self.mappings.append(best_mapping)
		    	
		    if not full_matched_transcripts:	
			# find events only for best transcript
			best_transcript = best_mapping.transcripts[0]
			events = self.find_events({best_transcript.id:all_block_matches[best_transcript.id]}, 
			                          align, 
			                          {best_transcript.id:best_transcript})
			for event in events:
			    event.contig_sizes.append(len(contig_seq))
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
		    fusion = fusion_finder.find_chimera(chimera_block_matches, transcripts, aligns, contig_seq, 
		                                        exon_bound_only=self.fusion_conditions['exon_bound_only'],
		                                        coding_only=self.fusion_conditions['coding_only'],
		                                        sense_only=self.fusion_conditions['sense_only'])
		    if fusion:
			homol_seq, homol_coords = None, None
			if self.aligner.lower() == 'gmap':
			    homol_seq, homol_coords = gmap.find_microhomology(alns[0], contig_seq)
			if homol_seq is not None:
			    fusion.homol_seq.append(homol_seq)
			    fusion.homol_coords.append(homol_coords)
			fusion.contig_sizes.append(len(contig_seq))
			self.events.append(fusion)
		
	# expand contig span
	for event in self.events:
	    # if contig_support_span is not defined (it can be pre-defined in ITD)
	    # then set it
	    if not event.contig_support_span:
		event.contig_support_span = event.contig_breaks
		if event.rearrangement == 'ins' or event.rearrangement == 'dup':
		    expanded_contig_breaks = expand_contig_breaks(event.chroms[0], 
			                                          event.breaks, 
			                                          event.contigs[0], 
			                                          [event.contig_breaks[0][0] + 1, event.contig_breaks[0][1] - 1], 
			                                          event.rearrangement, 
			                                          self.ref_fasta,
			                                          self.contigs_fasta,
			                                          self.debug)
		    if expanded_contig_breaks is not None:
			event.contig_support_span = [(expanded_contig_breaks[0] - 1, expanded_contig_breaks[1] + 1)]
					
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
		    		            
    def extract_novel_seq(self, adj):
	"""Extracts novel sequence in adjacency, should be a method of Adjacency"""
	aa = self.contigs_fasta.fetch(adj.contigs[0])
	contig_breaks = adj.contig_breaks[0]
	start, end = (contig_breaks[0], contig_breaks[1]) if contig_breaks[0] < contig_breaks[1] else (contig_breaks[1], contig_breaks[0])
	novel_seq = self.contigs_fasta.fetch(adj.contigs[0], start, end - 1)
	return novel_seq
	    
    def find_events(self, matches_by_transcript, align, transcripts, small=20):
	"""Find events from single alignment
	
	Wrapper for finding events within a single alignment
	Maybe a read-through fusion, calls fusion_finder.find_read_through()
	Maybe splicing or indels, call novel_splice_finder.find_novel_junctions()
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
	    small: (int) max size of ins or dup that the span of the novel seq would be used in contig_support_span,
	                 anything larger than that we just check if there is any spanning reads that cross the breakpoint
	Returns:
	    List of Adjacencies
	"""
	events = []
	
	# for detecting whether there is a read-through fusion
	genes = Set([transcripts[txt].gene for txt in matches_by_transcript.keys()])		
	if len(genes) > 1:
	    fusion = fusion_finder.find_read_through(matches_by_transcript, transcripts, align, 
	                                             exon_bound_only=self.fusion_conditions['exon_bound_only'],
	                                             coding_only=self.fusion_conditions['coding_only'],
	                                             sense_only=self.fusion_conditions['sense_only'])
	    if fusion is not None:
		events.append(fusion)
	    	    	
	# events within a gene
	local_events = novel_splice_finder.find_novel_junctions(matches_by_transcript, align, transcripts, self.ref_fasta)
	if local_events:
	    for event in local_events:
		adj = Event((align.target, align.target), event['pos'], '-', contig=align.query)
		if event['event'] in ('ins', 'del', 'dup', 'inv'):
		    adj.rearrangement = event['event']
		adj.rna_event = event['event']
		adj.genes = (list(genes)[0],)
		adj.transcripts = (event['transcript'][0],)
		adj.size = event['size']
		adj.orients = 'L', 'R'
		
		# converts exon index to exon number (which takes transcript strand into account)
		exons = event['exons'][0]
		adj.exons = map(transcripts[event['transcript'][0]].exon_num, exons)
		
		adj.contig_breaks = event['contig_breaks']
		
		if adj.rearrangement == 'ins' or adj.rearrangement == 'dup':
		    novel_seq = self.extract_novel_seq(adj)
		    adj.novel_seq = novel_seq if align.strand == '+' else reverse_complement(novel_seq)
		    
		    if len(novel_seq) >= self.itd_conditions['min_len']:
			itd_finder.detect_itd(adj, align, self.contigs_fasta.fetch(adj.contigs[0]), self.outdir, 
			                      self.itd_conditions['min_len'],
			                      self.itd_conditions['max_apart'],
			                      self.itd_conditions['min_pid'],
			                      debug=self.debug
			                      )
		    			
		    adj.size = len(adj.novel_seq)
		    
		    # bigger events mean we can only check if there is reads that cross the breakpoint
		    if adj.size > small:
			adj.contig_support_span = [(min(adj.contig_breaks[0]), min(adj.contig_breaks[0]) + 1)]
		    
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
    
    def find_support(self, r2c_bam_file=None, min_overlap=None, multi_mapped=False, perfect=True, get_seq=False, num_procs=1):
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
		#span = event.contig_breaks[i][0], event.contig_breaks[i][1]
		span = event.contig_support_span[i][0], event.contig_support_span[i][1]
		try:
		    coords[contig].append(span)
		except:
		    coords[contig] = [span]
		    
	support_reads = {}
	if coords:	    
	    avg_tlen = None
	    tlens = []
	    # if fewer than 10000 adjs, use 'fetch' of Pysam
	    if len(coords) < 1000 or num_procs == 1:
		support, tlens = fetch_support(coords, r2c_bam_file, self.contigs_fasta, 
		                               overlap_buffer=min_overlap, 
		                               perfect=perfect, 
		                               get_seq=get_seq, 
		                               debug=self.debug)
	    # otherwise use multi-process parsing of bam file
	    else:
		support, tlens = scan_all(coords, r2c_bam_file, self.contigs_fasta_file, num_procs, 
		                          overlap_buffer=min_overlap, 
		                          perfect=perfect, 
		                          get_seq=get_seq, 
		                          debug=self.debug)
		    
	    print 'total tlens', len(tlens)
		
	    for event in self.events:
		# initialize support
		event.support = {'spanning': 0, 'flanking': 0}
		support_contigs = []
		for i in range(len(event.contigs)):
		    contig = event.contigs[i]
		    span = event.contig_support_span[i][0], event.contig_support_span[i][1]
		    coord = '%s-%s' % (span[0], span[1])
					
		    if support.has_key(contig) and support[contig].has_key(coord):
			support_contigs.append(support[contig][coord][0])
			
			# if reads are to be extracted
			if get_seq:
			    if support[contig][coord][-1]:
				key = event.key(transcriptome=True)
				for read in support[contig][coord][-1]:
				    try:
					support_reads[key].add(read)
				    except:
					support_reads[key] = Set([read])
			
		if not multi_mapped:
		    event.support['spanning'] = sum(support_contigs)
		else:
		    event.support['spanning'] = max(support_contigs)
		    		     
	    if tlens:
		#avg_tlen = float(sum(tlens)) / len(tlens)
		print 'avg tlen', tlens[:10]
		    	
	if self.debug:
	    for event in self.events:
		print 'support', event.support
		
	return support_reads

    def screen_events(self, outdir, align_info=None, max_homol_allowed=None):
	"""Screen events identified and filter out bad ones
	
	Right now it just screens out fusion whose probe sequence can align to one single location
	
	Args:
	    outdir: (str) absolute path of output directory, for storing re-alignment results
	    aligner: (str) aligner name (gmap)
	    align_info: (dict) 'genome', 'index_dir', 'num_procs'
	"""
	bad_contigs = Set()
	if max_homol_allowed is not None:
	    for event in self.events:
		if event.homol_seq and len(event.homol_seq[0]) > max_homol_allowed:
		    if debug:
			sys.stdout.write('Screen out %s: homol_seq(%d-%d:%d) longer than maximum allowed(%d)\n' % (','.join(event.contigs),
			                                                                                     event.homol_coords[0][0],
			                                                                                     event.homol_coords[0][1],
			                                                                                     len(event.homol_seq[0]), 
			                                                                                     max_homol_allowed))
		    for contig in event.contigs:
			bad_contigs.add(contig)
	
	fusions = [e for e in self.events if e.rna_event == 'fusion']		
	if fusions and align_info is not None:
	    bad_contigs_realign = fusion_finder.screen_realigns(fusions, outdir, align_info, contigs_fasta=self.contigs_fasta, debug=self.debug)
	    if bad_contigs_realign:
		bad_contigs = bad_contigs.union(bad_contigs_realign)
		
	bad_event_indices = []
	if bad_contigs:
	    # remove any event that involve contigs that failed screening as mapping is not reliable
	    for e in reversed(range(len(self.events))):
		for contig in self.events[e].contigs:
		    if contig in bad_contigs and not e in bad_event_indices:
			bad_event_indices.append(e)
			break
		
	for e in bad_event_indices:
	    del self.events[e]
    
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
                    exon_bound_fusion_only=not options.include_non_exon_bound_fusion,
                    coding_fusion_only=not options.include_noncoding_fusion,
                    sense_fusion_only=not options.include_antisense_fusion,
                    debug=options.debug)
    em.map_contigs_to_transcripts()
	
    align_info = None
    if options.genome and options.index_dir and os.path.exists(options.index_dir):
	align_info = {
	    'aligner': em.aligner,
	    'genome': options.genome,
	    'index_dir': options.index_dir,
	    'num_procs': options.num_threads,
	}
    # screen events based on realignments
    em.screen_events(outdir, align_info=align_info, max_homol_allowed=options.max_homol_allowed)
    
    # added support
    if options.r2c_bam_file:
	support_reads = em.find_support(options.r2c_bam_file, options.min_overlap, options.multimapped, 
	                                num_procs=options.num_threads, 
	                                get_seq=options.output_support_reads)
	
    # merge events captured by different contigs into single events
    events_merged = Adjacency.merge(em.events, transcriptome=True)
    
    # filter read support after merging
    if options.r2c_bam_file:
	Event.filter_by_support(events_merged, options.min_support)
	
    # output events
    Event.output(events_merged, outdir, sort_by_event_type=options.sort_by_event_type)
    
    # output support reads
    if options.output_support_reads and support_reads:
	Event.output_reads(events_merged, support_reads, '%s/support_reads.fa' % outdir)
    
    # output mappings
    Mapping.output(em.mappings, '%s/contig_mappings.tsv' % outdir)
    gene_mappings = Mapping.group(em.mappings)
    Mapping.output(gene_mappings, '%s/gene_mappings.tsv' % outdir)
    
if __name__ == '__main__':
    usage = "Usage: %prog c2g_bam aligner contigs_fasta annotation_file genome_file(indexed) out_dir"
    parser = OptionParser(usage=usage, version=__version__)
    
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
    parser.add_option("--min_support", dest="min_support", help="minimum read support. Default:2", type='int', default=2)
    parser.add_option("--include_non_exon_bound_fusion", dest="include_non_exon_bound_fusion", help="include fusions where breakpoints are not at exon boundaries", 
                      action="store_true", default=False)
    parser.add_option("--include_noncoding_fusion", dest="include_noncoding_fusion", help="include non-coding genes in detecting fusions", action="store_true", default=False)
    parser.add_option("--include_antisense_fusion", dest="include_antisense_fusion", help="include antisense fusions", action="store_true", default=False)
    parser.add_option("--sort_by_event_type", dest="sort_by_event_type", help="sort output by event type", action="store_true", default=False)
    parser.add_option("--output_support_reads", dest="output_support_reads", help="output support reads", action="store_true", default=False)
    parser.add_option("--max_homol_allowed", dest="max_homol_allowed", help="maximun amount of microhomology allowed. Default:10", type="int", default=10)
    parser.add_option("--debug", dest="debug", help="debug mode", action="store_true", default=False)
    
    (options, args) = parser.parse_args()
    if len(args) == 6:
        main(args, options)     