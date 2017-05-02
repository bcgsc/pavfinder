import argparse
import sys
import os
import pysam
from intspan import intspan
from itertools import groupby
from sets import Set
from alignment import Alignment
from transcript import Transcript
from adjacency import Adjacency
from collections import defaultdict
from read_support import find_support
from novel_splice_finder import find_novel_junctions, report_items, report, extract_features, filter_events

class ExonMapper:
    def __init__(self, annot, transcripts_dict, genome_fasta,
                 suppl_annots = None,
                 debug = False):
	self.transcripts_dict = transcripts_dict
	self.genome_fasta = genome_fasta
	self.annot = annot
	self.debug = debug
        
	self.suppl_annots = suppl_annots

    def map_aligns(self, bam, query_fasta, genome_fasta, accessory_known_features=None, find_events=True):
	mappings = defaultdict(list)
	junc_adjs = []
	events = []
	for query, group in groupby(bam.fetch(until_eof=True), lambda aln: aln.query_name):
	    print 'processing', query
	    aligns = []
	    for aln in list(group):
		if not aln.is_unmapped:
		    aligns.append(Alignment.from_alignedRead(aln, bam))
		
	    if not aligns:
		continue
	    
	    query_seq = query_fasta.fetch(query)
	    	    
	    for align in aligns:
		if not align.has_canonical_target():
		    continue
		block_matches = self.map_align(align)
		if block_matches:
		    tid = self.pick_best_mapping(block_matches, align)
		    if tid is not None:
			transcript = self.transcripts_dict[tid]
			olap = self.overlap(align, transcript)
			mappings[query].append((transcript.gene, transcript.id, olap))
			
			junc_adjs.extend(self.collect_junctions(align, transcript, block_matches[tid]))
			
			if find_events:
			    events.extend(find_novel_junctions(block_matches[tid],
			                                       align,
			                                       transcript,
			                                       query_seq,
			                                       self.genome_fasta,
			                                       accessory_known_features=accessory_known_features)
			                  )
	
	return mappings, junc_adjs, events

    @classmethod
    def sort_adjs(cls, adjs):
	dict_for_sorting = {}
	for i in range(len(adjs)):
	    key = adjs[i].chroms[0], adjs[i].genome_breaks[0], adjs[i].genome_breaks[1]
	    dict_for_sorting[key] = i

	adjs_sorted = []
	for key in sorted(dict_for_sorting.keys(), Adjacency.cmp_genome_coords):
	    adjs_sorted.append(adjs[dict_for_sorting[key]])

	return adjs_sorted

    @classmethod
    def collect_junctions(cls, align, transcript, block_matches):
	adjs = []
	for i in range(len(align.blocks) - 1):
	    if block_matches[i] is None or block_matches[i+1] is None:
		continue
	    
	    if block_matches[i][-1][0] == block_matches[i + 1][0][0]:
		continue
	    
	    adj = Adjacency(align.query,
	                    (align.target, align.target),
	                    (align.query_blocks[i][1], align.query_blocks[i + 1][0]),
	                    (align.blocks[i][1], align.blocks[i + 1][0])
	                    )
	    adj.genome_breaks = adj.target_breaks
	    adj.chroms = adj.targets
	    adj.transcripts = (transcript, transcript)
	    adj.exons = (transcript.exon_num(block_matches[i][-1][0]),
	                 transcript.exon_num(block_matches[i + 1][0][0]))
	    adjs.append(adj)

	return adjs
    
    @classmethod
    def output_mappings(cls, mappings, out_file):
	# group by gene ...
	groups_by_gene = defaultdict(list)
	for query in sorted(mappings.keys()):
	    for mapping in mappings[query]:
		groups_by_gene[mapping[0]].append((query, mapping[0], mapping[1], mapping[2]))

	# ... and sort by coverage
	out = open(out_file, 'w')
	out.write('%s\n' % '\t'.join(('contig', 'gene', 'transcript', 'coverage')))
	for gene in sorted(groups_by_gene.keys()):
	    for mapping in sorted(groups_by_gene[gene], key = lambda m: m[3], reverse=True):
		query, gene, tid, olap = mapping
		out.write('%s\n' % '\t'.join((query, gene, tid, '%.3f' % olap)))
	out.close()
	
    @classmethod
    def output_juncs(cls, juncs, out_file):
	def adj_to_bed(adj):
	    cols = []
	    cols.append(adj.chroms[0])
	    cols.append(adj.genome_breaks[0] - 1)
	    cols.append(adj.genome_breaks[1])
	    label = [adj.transcripts[0].gene,
	             adj.transcripts[0].id,
	             'E%d' % adj.exons[0],
	             'E%d' % adj.exons[1],
	             ]
	    cols.append('.'.join(map(str, label)))

	    if adj.spanning is not None:
		cols.append(adj.spanning)
	    else:
		cols.append('.')
	    if adj.transcripts[0].strand == '+':
		cols.append('+')
	    else:
		cols.append('-')

	    cols.append(adj.genome_breaks[0] - 1)
	    cols.append(adj.genome_breaks[1])
	    cols.append(0)
	    cols.append(2)
	    cols.append('1,1')
	    cols.append('0,%s' % (adj.genome_breaks[1] - adj.genome_breaks[0]))
	    return '\t'.join(map(str, cols))
	    
	out = open(out_file, 'w')
	for junc in cls.sort_adjs(juncs):
	    out.write('%s\n' % adj_to_bed(junc))
	out.close()
	
    @classmethod
    def output_events(cls, events, out_file, header=None):
	def create_bedpe_header():
	    cols = []
	    for i in {1,2}:
		for label in {'chrom', 'start', 'end'}:
		    cols.append('%s%d' % (label, i))
	    for label in ('name', 'score', 'strand1', 'strand2'):
		cols.append(label)

	    return cols + report_items.keys()

	def group_events():
	    grouped_events = defaultdict(list)
	    for event in events:
		grouped_events[event.event].append(event)
		
	    events_ordered = []
	    for event_type in event_types:
		sorted_events = cls.sort_adjs(grouped_events[event_type])
		if event_type in ('retained_intron', 'novel_exon'):
		    indices_ordered = []
		    for i in range(len(sorted_events)):
			if i in indices_ordered:
			    continue
			event = sorted_events[i]
			
			for j in range(i + 1, len(sorted_events)):
			    if j in indices_ordered:
				continue
			    if sorted_events[j].link == event:
				indices_ordered.append(i)
				indices_ordered.append(j)
		    for i in indices_ordered:
			events_ordered.append(sorted_events[i])
		else:
		    for event in sorted_events:
			events_ordered.append(event)

	    return events_ordered
 
	event_types = ['skipped_exon',
	               'novel_exon',
	               'novel_intron',
	               'novel_acceptor',
	               'novel_donor',
	               'retained_intron'
	               ]
	out = open(out_file, 'w')
	if header is not None:
	    if type(header) is str:
		out.write('#%s\n' % header)
	    elif type(header) is tuple or type(header) is list:
		for h in header:
		    out.write('#%s\n' % h)
	out.write('%s\n' % '\t'.join(create_bedpe_header()))
	events_ordered = group_events()

	counter = 1
	for i in range(len(events_ordered)):
	    if events_ordered[i].event in ('retained_intron', 'novel_exon'):
		if i > 0 and\
		   events_ordered[i - 1].event == events_ordered[i].event and\
		   events_ordered[i - 1].seq_id == events_ordered[i].seq_id:
		    event_id = '%d.2' % counter
		    counter += 1
		else:
		    event_id = '%d.1' % counter
	    else:
		event_id = counter

	    out.write('%s\n' % report(events_ordered[i], event_id=event_id))
	    
	    if events_ordered[i].event not in ('retained_intron', 'novel_exon'):
		counter += 1
	out.close()
	
    def overlap(self, align, transcript):
        """Overlaps alignment spans with exon-spans of each matching transcript"""
	def create_span(blocks):
	    """Creates intspan for each block"""
	    span = None
	    for block in blocks:
		if (type(block) is tuple or type(block) is list) and len(block) == 2:
		    if span is None:
			span = intspan('%s-%s' % (block[0], block[1]))
		    else:
			span = span.union(intspan('%s-%s' % (block[0], block[1])))
	    return span

        align_span = create_span(align.blocks)
	exon_span = create_span(transcript.exons)
	olap = exon_span.intersection(align_span)
	return float(len(olap)) / float(len(exon_span))
	    
    def map_align(self, align):
	mappings = {}
	for record in self.annot.fetch(align.target, align.tstart, align.tend):
	    if not self.transcripts_dict.has_key(record.transcript_id):
		continue

	    transcript = self.transcripts_dict[record.transcript_id]	    
	    mapping = self.map_exons(align.blocks, transcript.exons)
	    if len([m for m in mapping if m is None]) != len(mapping):
		mappings[transcript.id] = mapping
	    
	return mappings
	    
    def map_exons(self, blocks, exons):
	"""Maps alignment blocks to exons
		
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
	mappings = []
	for b in range(len(blocks)):
	    block_matches = []
	    
	    for e in range(len(exons)):
		block_match = self.match_exon(blocks[b], exons[e])
		if block_match != '':
		    block_matches.append((e, block_match))
		    
	    if not block_matches:
		block_matches = None
	    mappings.append(block_matches)
		
	return mappings
    
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
        match = ''
        
	if min(block[1], exon[1]) - max(block[0], exon[0]) > 0:
	    for i in range(0, 2):
		if block[i] == exon[i]:
		    match += '='
		elif block[i] > exon[i]:
		    match += '>'
		else:
		    match += '<'
        
        return match
    
    def pick_best_mapping(self, mappings, align, debug=False):
	"""Selects best mapping among transcripts"""
	def calc_score(matches):
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

	    return score
	
	def calc_penalty(matches):
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
			penalty += 1

	    return penalty
	
	def calc_from_edge(matches, transcript):
	    from_edge = None
	    if matches[0] is None or matches[-1] is None:
		from_edge = 10000
	    else:
		if transcript.strand == '+':
		    start_exon_num = matches[0][0][0] + 1
		    end_exon_num = matches[-1][-1][0] + 1
		else:
		    start_exon_num = transcript.num_exons() - matches[0][0][0]
		    end_exon_num = transcript.num_exons() - matches[-1][-1][0]
		    
		from_edge = align.tstart - transcript.exon(start_exon_num)[0] + \
		          align.tend - transcript.exon(end_exon_num)[1]
	    return from_edge

	metrics = {}
	for tid, matches in mappings.iteritems():
	    transcript = self.transcripts_dict[tid]
	    metric = {}
	    # points are scored for matching exon boundaries
	    metric['score'] = calc_score(matches) - calc_penalty(matches)
	    metric['from_edge'] = calc_from_edge(matches, self.transcripts_dict[tid])
	    metric['txt_size'] = transcript.length()
	    metrics[tid] = metric
	    	    
	    if debug:
		sys.stdout.write("mapping %s %s %s\n" % (align.query, tid, metric))
	    
	if mappings:
	    tids_sorted = sorted(metrics.keys(), key = lambda tid: (-1 * metrics[tid]['score'], 
		                                                    metrics[tid]['from_edge'], 
		                                                    metrics[tid]['txt_size']))
	    if debug:
		for tid in tids_sorted:
		    sys.stdout.write('sorted %s %s\n' % (tid, metrics[tid]))
		    
	    return tids_sorted[0]
	return None
    
    def is_full_match(self, block_matches, perfect=True):
	"""Determines if the query 'fully' matches

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
	if len(block_matches) == 1:
	    if block_matches[0] is not None and len(block_matches[0]) == 1:
		if block_matches[0][0][1] == '==' or block_matches[0][0][1] == '>=' or block_matches[0][0][1] == '=<':
		    return True
	    
	    return False
	
	# if a block is mapped to >1 exon
	if None in block_matches or [m for m in block_matches if len(m) > 1]:
	    return False
	
	exons = [m[0][0] for m in block_matches]
	
	# the last condition is to check if the exons are consecutive
	if block_matches[0][0][1][1] == '=' and\
	   block_matches[-1][0][1][0] == '=' and\
	   len([m for m in block_matches[1:-1] if m[0][1] == '==']) == len(block_matches) - 2 and\
	   (len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(block_matches) - 1 or\
	    len([(a, b) for a, b in zip(exons, exons[1:]) if b == a - 1]) == len(block_matches) - 1):
	    return True
	if block_matches[0][0][1][1] == '=' and\
	   block_matches[-1][0][1][0] == '=' and\
	   len([m for m in block_matches[1:-1] if m[0][1] == '==']) == len(block_matches) - 2:
	    if perfect:
		if (len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(block_matches) - 1 or\
		    len([(a, b) for a, b in zip(exons, exons[1:]) if b == a - 1]) == len(block_matches) - 1):
		    return True
	    else:
		return True
	
	return False
