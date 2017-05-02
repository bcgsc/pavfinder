from sets import Set
import sys
import re
from collections import defaultdict, OrderedDict
from alignment import reverse_complement
from copy import deepcopy

class Adjacency:

    report_items = OrderedDict(
        [('ID', None),
         ('event', 'event'),
         ('chrom1', 'chroms,0'),
         ('genome_break1', 'genome_breaks,0'),
         ('orient1', 'orients,0'),
         ('chrom2', 'chroms,1'),
         ('genome_break2', 'genome_breaks,1'),
         ('orient2', 'orients,1'),
         ('size', 'size'),
         ('gene1', None),
         ('transcript1', None),
         ('transcript_break1', 'transcript_breaks,0'),
         ('exon1', 'exons,0',),
         ('exon_bound1', 'exon_bounds,0'),
         ('gene2', None),
         ('transcript2', None),
         ('transcript_break2', 'transcript_breaks,1'),
         ('exon2', 'exons,1',),
         ('exon_bound2', 'exon_bounds,1'),
         ('gene_5prime', None),
         ('gene_3prime', None),
         ('exon_5prime', 'exons_oriented,0'),
         ('exon_3prime', 'exons_oriented,1'),
         ('feature', 'feature'),
         ('seq_id', 'seq_id'),
         ('seq_breaks', 'seq_breaks'),
         ('ins_seq', 'ins_seq'),
         ('homol_seq', 'homol_seq'),
         ('homol_seq_coords', 'homol_seq_coords'),
         ('novel_seq', 'novel_seq'),
         ('novel_seq_coords', 'novel_seq_coords'),
         ('copy_number_change', None),
         ('repeat_seq', 'repeat_seq'),
         ('in_frame', 'in_frame'),
         ('probe', 'probe'),
         ('support_span', 'support_span'),
         ('spanning_reads', 'spanning'),
         ('flanking_pairs', 'flanking'),
         ]
    )
    
    bedpe_items = OrderedDict(
        [('chrom1', 'chroms,0'),
         ('start1', None),
         ('end1', 'genome_breaks,0'),
         ('chrom2', 'chroms,1'),
         ('start2', None),
         ('end2', 'genome_breaks,1'),
         ('name', None),
         ('score', None),
         ('strand1', None),
         ('strand2', None),
         ('orient1', 'orients,0'),
         ('orient2', 'orients,1'),
         ('event', 'event'),
         ('size', 'size'),
         ('gene1', None),
         ('transcript1', None),
         ('transcript_break1', 'transcript_breaks,0'),
         ('exon1', 'exons,0',),
         ('exon_bound1', 'exon_bounds,0'),
         ('gene2', None),
         ('transcript2', None),
         ('transcript_break2', 'transcript_breaks,1'),
         ('exon2', 'exons,1',),
         ('exon_bound2', 'exon_bounds,1'),
         ('gene_5prime', None),
         ('gene_3prime', None),
         ('exon_5prime', 'exons_oriented,0'),
         ('exon_3prime', 'exons_oriented,1'),
         ('feature', 'feature'),
         ('seq_id', 'seq_id'),
         ('seq_breaks', 'seq_breaks'),
         ('ins_seq', 'ins_seq'),
         ('homol_seq', 'homol_seq'),
         ('homol_seq_coords', 'homol_seq_coords'),
         ('novel_seq', 'novel_seq'),
         ('novel_seq_coords', 'novel_seq_coords'),
         ('copy_number_change', None),
         ('repeat_seq', 'repeat_seq'),
         ('in_frame', 'in_frame'),
         ('probe', 'probe'),
         ('support_span', 'support_span'),
         ('spanning_reads', 'spanning'),
         ('flanking_pairs', 'flanking'),
         ]
    )

    event_types = ['fusion', 'read_through', 'ITD', 'PTD', 'dup', 'ins', 'del', 'dup_inv', 'repeat_expansion', 'repeat_reduction']

    def __init__(self, seq_id, targets, seq_breaks, target_breaks,
                 chroms = None, genome_breaks = None, orients = None,
                 transcripts = None, transcript_breaks = None,
                 exons = None, exon_bounds = None,
                 rearrangement = None, event = None,
                 novel_seq = None, novel_seq_coords = None,
                 homol_seq = None, homol_seq_coords = None,
                 ins_seq = None, ins_seq_coords = None,
                 copy_num_change = None, repeat_seq = None,
                 upstream_transcript = None, downstream_transcript = None,
                 exons_oriented = None, exon_bounds_oriented = None,
                 feature = None, effect = None, in_frame = None, sense_fusion=None,
                 support_span=None, spanning=None, flanking=None,
                 probe=None, size=None, link=None):
	for attr, value in locals().iteritems():
	    setattr(self, attr, value)
    
    def update_attrs(self, **kwargs):
	for attr, value in kwargs.iteritems():
	    if hasattr(self, attr):
		existing_value = getattr(self, attr)
		if existing_value is None:
		    setattr(self, attr, value)
		
    def update_genome_breaks(self):
	if self.transcript_breaks and self.transcripts and not self.genome_breaks:
	    genome_breaks = []
	    for i in range(len(self.transcript_breaks)):
		genome_breaks.append(self.transcripts[i].txt_coord_to_genome_coord(self.transcript_breaks[i]))
	    self.genome_breaks = tuple(genome_breaks)
	    
    def update_transcript_breaks(self):
	if self.genome_breaks and self.transcripts and not self.transcript_breaks:
	    transcript_breaks = []
	    for i in range(len(self.genome_breaks)):
		transcript_breaks.append(self.transcripts[i].genome_coord_to_txt_coord(self.genome_breaks[i]))
	    self.transcript_breaks = tuple(transcript_breaks)
	    
    def update_exons(self, target_type):
	"""Update exons and exon_bounds

	Argument:
	    target_type: "transcripts" or "genome"
	    if target_type is "genome", orient and transcript.strand need to be considered
	    if target_type is "transcripts", only orient needs to be considered
	"""
	def check_exon_bound(exon_bound_edge, orient, transcript):
	    exon_bound = False
	    if exon_bound_edge is not None:
		if target_type == 'genome':
		    if transcript.strand == '+':
			if (orient == 'L' and exon_bound_edge == 3) or\
			   (orient == 'R' and exon_bound_edge == 5):
			    exon_bound = True
		    else:
			if (orient == 'R' and exon_bound_edge == 3) or\
			   (orient == 'L' and exon_bound_edge == 5):
			    exon_bound = True
		elif target_type == 'transcripts':
		    if (orient == 'L' and exon_bound_edge == 3) or\
		       (orient == 'R' and exon_bound_edge == 5):
			exon_bound = True
	    return exon_bound

	if not self.exons and not self.exon_bounds:
	    exons = []
	    exon_bounds = []
	    if self.transcripts:
		if self.transcript_breaks:
		    for i in range(len(self.transcript_breaks)):
			exons.append(self.transcripts[i].txt_coord_to_exon(self.transcript_breaks[i]))
			exon_bound_edge = self.transcripts[i].at_exon_bound(txt_coord = self.transcript_breaks[i],
			                                                    exon_num = exons[-1],
			                                                    return_edge=True)
			exon_bounds.append(check_exon_bound(exon_bound_edge, self.orients[i], self.transcripts[i]))

		elif self.genome_breaks:
		    for i in range(len(self.genome_breaks)):
			exons.append(self.transcripts[i].coord_to_exon(self.genome_breaks[i]))
			exon_bound_edge = self.transcripts[i].at_exon_bound(genome_coord = self.genome_breaks[i],
			                                                    exon_num = exons[-1],
			                                                    return_edge = True)
			exon_bounds.append(check_exon_bound(exon_bound_edge, self.orients[i], self.transcripts[i]))

		if exons and exon_bounds:
		    self.exons = tuple(exons)
		    self.exon_bounds = tuple(exon_bounds)

    def change_tuple(self, attr, index, value):
	"""Change value in tuple-type attribute (used in adjust_transcript_breaks()"""
	values = list(getattr(self, attr))
	values[index] = value
	setattr(self, attr, tuple(values))

    def adjust_transcript_breaks(self, align_strands):
	"""Adjust breakpoint positions if exon_bound not True when there is homology"""
	def check_exon_bound(exon_bound_edge, orient):
	    exon_bound = False
	    if exon_bound_edge is not None:
		if (orient == 'L' and exon_bound_edge == 3) or\
		   (orient == 'R' and exon_bound_edge == 5):
		    exon_bound = True
	    return exon_bound

	seq_breaks = self.seq_breaks
	transcript_breaks = self.transcript_breaks

	if not self.exon_bounds[0]:
	    for i in range(1, len(self.homol_seq) + 1):
		sbreak = seq_breaks[0] - i
		if align_strands[0] == '+':
		    tbreak = transcript_breaks[0] - i
		else:
		    tbreak = transcript_breaks[0] + i
		homol_seq = self.homol_seq[:-1*i]
		homol_seq_coords = list(self.homol_seq_coords)
		homol_seq_coords[1] -= i
		if not homol_seq:
		    homol_seq = None
		    homol_seq_coords = None
		exon_num = self.transcripts[0].txt_coord_to_exon(tbreak)
		exon_bound_edge = self.transcripts[0].at_exon_bound(txt_coord = tbreak,
		                                                    exon_num = exon_num,
		                                                    return_edge=True)
		gbreak = self.transcripts[0].txt_coord_to_genome_coord(tbreak)

		if check_exon_bound(exon_bound_edge, self.orients[0]):
		    self.change_tuple('exons', 0, exon_num)
		    self.change_tuple('exon_bounds', 0, True)
		    self.change_tuple('transcript_breaks', 0, tbreak)
		    self.change_tuple('genome_breaks', 0, gbreak)
		    self.change_tuple('seq_breaks', 0, sbreak)
		    self.homol_seq = homol_seq
		    if type(homol_seq_coords) is list:
			self.homol_seq_coords = tuple(homol_seq_coords)
		    else:
			self.homol_seq_coords = homol_seq_coords
		    break

	if self.homol_seq and not self.exon_bounds[1]:
	    for i in range(1, len(self.homol_seq) + 1):
		sbreak = seq_breaks[1] + i
		if align_strands[1] == '+':
		    tbreak = transcript_breaks[1] + i
		else:
		    tbreak = transcript_breaks[1] - i
		homol_seq = self.homol_seq[i:]
		homol_seq_coords = list(self.homol_seq_coords)
		homol_seq_coords[0] += i
		if not homol_seq:
		    homol_seq = None
		    homol_seq_coords = None
		exon_num = self.transcripts[1].txt_coord_to_exon(tbreak)
		exon_bound_edge = self.transcripts[1].at_exon_bound(txt_coord = tbreak,
		                                                    exon_num = exon_num,
		                                                    return_edge=True)
		gbreak = self.transcripts[1].txt_coord_to_genome_coord(tbreak)

		if check_exon_bound(exon_bound_edge, self.orients[1]):
		    self.change_tuple('exons', 1, exon_num)
		    self.change_tuple('exon_bounds', 1, True)
		    self.change_tuple('transcript_breaks', 1, tbreak)
		    self.change_tuple('genome_breaks', 1, gbreak)
		    self.change_tuple('seq_breaks', 1, sbreak)
		    self.homol_seq = homol_seq
		    if type(homol_seq_coords) is list:
			self.homol_seq_coords = tuple(homol_seq_coords)
		    else:
			self.homol_seq_coords = homol_seq_coords
		    break

    def details(self):
	attr_values = []
	for attr, value in self.__dict__.iteritems():
	    attr_values.append('%s:%s' % (attr, value))
	    
	transcript_ids = ['NA', 'NA']
	genes = ['NA', 'NA']
	genes_ordered = ['NA', 'NA']
	if self.transcripts and len(self.transcripts) == 2:
	    transcript_ids = (self.transcripts[0].id, self.transcripts[1].id)
	    genes = (self.transcripts[0].gene, self.transcripts[1].gene)
	if self.upstream_transcript and self.downstream_transcript:
	    genes_ordered = (self.upstream_transcript.gene, self.downstream_transcript.gene)
	return ' '.join(attr_values) +\
	       ' tids:%s' % ','.join(transcript_ids) +\
	       ' genes:%s' % ','.join(genes) +\
	       ' genes_ordered:%s' % ','.join(genes_ordered)
	
    def get_size(self):
	size = 'na'
	
	homol_seq_len = 0
	if self.homol_seq is not None and self.homol_seq != '-':
	    homol_seq_len = len(self.homol_seq)

	novel_seq_len = 0
	if self.novel_seq is not None and self.novel_seq != '-':
	    novel_seq_len = len(self.novel_seq)
	    
	if not self.transcript_breaks or not self.seq_breaks or\
	   self.transcripts[0] != self.transcripts[1]:
	    return size
	    
	if self.rearrangement == 'dup' or self.rearrangement == 'inv-dup':
	    if self.transcript_breaks[0] is not None and self.transcript_breaks[1] is not None:
		breaks_sorted = sorted(self.transcript_breaks)
		size = breaks_sorted[1] - breaks_sorted[0] + 1 + novel_seq_len - homol_seq_len
	    
	elif self.rearrangement == 'ins':
	    if self.seq_breaks[0] is not None and self.seq_breaks[1] is not None:
		seq_breaks_sorted = sorted(self.seq_breaks)
		size = seq_breaks_sorted[1] - seq_breaks_sorted[0] - 1 + novel_seq_len - homol_seq_len
	    
	elif self.rearrangement == 'del':
	    if self.target_breaks[0] is not None and self.target_breaks[1] is not None:
		target_breaks_sorted = sorted(self.target_breaks)
		size = target_breaks_sorted[1] - target_breaks_sorted[0] - 1 + homol_seq_len
        
        return size
    
    def update_support_span(self):
	sorted_seq_breaks = sorted(self.seq_breaks)
	if self.support_span is None:
	    self.support_span = sorted_seq_breaks
	    
    def report(self, event_id=None):
	data = []
	for item, label in self.report_items.iteritems():
	    value = 'na'
	    if item == "name":
		if event_id is not None:
		    value = event_id
		    
	    elif label is not None and ',' in label:
		attr, index = label.split(',')
		if hasattr(self, attr):
		    values = getattr(self, attr)
		    if (type(values) is tuple or type(values) is list) and\
		       len(values) > int(index):
			value = values[int(index)]
			    
	    elif item[:4] == 'gene' and item[-1].isdigit():
		if item[-1] == '1':
		    value = self.transcripts[0].gene
		elif item[-1] == '2':
		    value = self.transcripts[1].gene
		    
	    elif item == 'gene_5prime':
		if self.upstream_transcript is not None:
		    value = self.upstream_transcript.gene
		else:
		    value = 'na'
		
	    elif item == 'gene_3prime':
		if self.downstream_transcript is not None:
		    value = self.downstream_transcript.gene
		else:
		    value = 'na'
		
	    elif item[:len('transcript')] == 'transcript' and item[-1].isdigit():
		if item[-1] == '1':
		    value = self.transcripts[0].id
		elif item[-1] == '2':
		    value = self.transcripts[1].id
		    
	    elif item == 'copy_number_change':
		if self.copy_num_change is not None:
		    value = '>'.join(map(str, self.copy_num_change))
		else:
		    value = 'na'
			
	    elif hasattr(self, label):
		val = getattr(self, label)
		if val is not None:
		    if (type(val) is tuple or type(val) is list) and\
		       len(val) == 2:
			value = '%s-%s' % (val[0], val[1])
		    else:
			value = val
		    
	    data.append(str(value))

	return '\t'.join(data)

    def as_bedpe(self, event_id=None):
	data = []
	for item, label in self.bedpe_items.iteritems():
	    value = 'na'
	    if item == "name":
		if event_id is not None:
		    value = event_id
		else:
		    value = '.'

	    elif item == 'start1':
		value = self.genome_breaks[0] - 1

	    elif item == 'start2':
		value = self.genome_breaks[1] - 1

	    elif 'strand' in item:
		value = '+'

	    elif item == 'score':
		value = '.'

	    elif label is not None and ',' in label:
		attr, index = label.split(',')
		if hasattr(self, attr):
		    values = getattr(self, attr)
		    if (type(values) is tuple or type(values) is list) and\
		       len(values) > int(index):
			value = values[int(index)]

	    elif item[:4] == 'gene' and item[-1].isdigit():
		if item[-1] == '1':
		    value = self.transcripts[0].gene
		elif item[-1] == '2':
		    value = self.transcripts[1].gene

	    elif item == 'gene_5prime':
		if self.upstream_transcript is not None:
		    value = self.upstream_transcript.gene
		else:
		    value = 'na'

	    elif item == 'gene_3prime':
		if self.downstream_transcript is not None:
		    value = self.downstream_transcript.gene
		else:
		    value = 'na'

	    elif item[:len('transcript')] == 'transcript' and item[-1].isdigit():
		if item[-1] == '1':
		    value = self.transcripts[0].id
		elif item[-1] == '2':
		    value = self.transcripts[1].id

	    elif item == 'copy_number_change':
		if self.copy_num_change is not None:
		    value = '>'.join(map(str, self.copy_num_change))
		else:
		    value = 'na'

	    elif hasattr(self, label):
		val = getattr(self, label)
		if val is not None:
		    if (type(val) is tuple or type(val) is list) and\
		       len(val) == 2:
			value = '%s-%s' % (val[0], val[1])
		    else:
			value = val

	    data.append(str(value))

	return '\t'.join(data)

    @classmethod
    def report_events(cls, events, outfile, sort_by_coord=False, header=None):
	def compare_event(e1, e2):
	    cmp_coord_result = cls.cmp_genome_coords((e1.chroms[0], e1.genome_breaks[0], e1.genome_breaks[1]),
	                                             (e2.chroms[0], e2.genome_breaks[0], e2.genome_breaks[1]))
	    if cmp_coord_result != 0:
		return cmp_coord_result

	    cmp_coord_result = cls.cmp_genome_coords((e1.chroms[1], e1.genome_breaks[1], e1.genome_breaks[1]),
	                                             (e2.chroms[1], e2.genome_breaks[1], e2.genome_breaks[1]))

	    if cmp_coord_result != 0:
		return cmp_coord_result

	    if e1.support is not None and e2.support is not None:
		if int(e1.support) > int(e2.support):
		    return -1
		elif int(e2.support) > int(e1.support):
		    return 1

	    return 0

	out = open(outfile, 'w')
	if header is not None:
	    if type(header) is str:
		out.write('#%s\n' % header)
	    elif type(header) is tuple or type(header) is list:
		for h in header:
		    out.write('#%s\n' % h)
	out.write('#%s\n' % '\t'.join(cls.bedpe_items.keys()))
	
	if not sort_by_coord:
	    events_grouped = defaultdict(list)
	    for event in events:
		events_grouped[event.event].append(event)
	    others = []
	    counter = 1
	    for event_type in cls.event_types:
		if events_grouped.has_key(event_type):
		    events_sorted = sorted(events_grouped[event_type], cmp = lambda e1, e2 : compare_event(e1, e2))
		    for event in events_sorted:
			out.write('%s\n' % event.as_bedpe(event_id=counter))
			counter += 1
		else:
		    others.append(event_type)
	    for event_type in others:
		events_sorted = sorted(events_grouped[event_type], cmp = lambda e1, e2 : compare_event(e1, e2))
		for event in events_sorted:
		    out.write('%s\n' % event.as_bedpe(event_id=counter))
		    counter += 1

	else:
	    events_sorted = sorted(events, cmp = lambda e1, e2 : compare_event(e1, e2))
	    counter = 1
	    for event in events_sorted:
		out.write('%s\n' % event.as_bedpe(event_id=counter))
		counter += 1

        out.close()
    
    def as_tab(self):        
        outputs = []
	
	data = []
	data.append(self.id)
	data.append(self.rearrangement)
	for j in range(2):
	    data.append(self.chroms[j])
	    data.append(str(self.breaks[j]))
	    data.append(self.orients[j])
	data.append(str(self.get_size()))
	data.append(','.join(self.contigs))
	data.append(','.join([str(b[0]) for b in self.contig_breaks]))
	data.append(','.join([str(b[1]) for b in self.contig_breaks]))
		
	if self.homol_seq:
	    data.append(','.join(self.homol_seq))
	else:
	    data.append('-')
	    
	if self.homol_coords:
	    data.append(','.join([str(b[0]) for b in self.homol_coords]))
	    data.append(','.join([str(b[1]) for b in self.homol_coords]))
	else:
	    data.append('-')
	    data.append('-')
	    
	data.append(self.novel_seq)
	
	try:
	    data.append(self.probes[0])
	except:
	    data.append('-')
	if self.support is not None:
	    data.append(str(self.support['spanning']))
	else:
	    data.append('-')
	outputs.append('\t'.join(map(str, data)))
	                
        return '\n'.join(outputs)
        
    def get_contig_support_span(self, contig_index):
	try:
	    return (self.homol_coords[contig_index][0], self.homol_coords[contig_index][1])
	except:
	    return (self.contig_breaks[contig_index][0], self.contig_breaks[contig_index][1])

    @classmethod
    def cmp_genome_coords(cls, coord1, coord2):
	def extract_chrom(chrom):
	    if chrom.isdigit():
		return None, int(chrom)
	    else:
		m = re.search(r'^(\D+)(\d+)$', chrom)
		if m:
		    return m.group(1), int(m.group(2))
	    return None, chrom

	if coord1[0] != coord2[0]:
	    prefix1, chrom1 = extract_chrom(coord1[0])
	    prefix2, chrom2 = extract_chrom(coord2[0])

	    if prefix1 == prefix2:
		if chrom1 < chrom2:
		    return -1
		else:
		    return 1
	    else:
		if coord1[0] < coord2[0]:
		    return -1
		else:
		    return 1
	else:
	    if int(coord1[1]) < int(coord2[1]):
		return -1
	    elif int(coord1[1]) > int(coord2[1]):
		return 1
	    elif int(coord1[2]) < int(coord2[2]):
		return -1
	    elif int(coord1[2]) > int(coord2[2]):
		return 1
	    else:
		return 0

    def is_genome_breaks_sorted(self):
	if self.cmp_genome_coords((self.chroms[0], self.genome_breaks[0], self.genome_breaks[1]),
	                          (self.chroms[1], self.genome_breaks[1], self.genome_breaks[1])) <= 0:
	    return True
	else:
	    return False

    def reverse_genome_breaks(self, change=False, annotation=False):
	def flip(fields, new_vals):
	    for field in fields:
		new_vals[field] = None
		if hasattr(self, field):
		    val = getattr(self, field)
		    if type(val) is tuple or type(val) is list:
			new_vals[field] = tuple(reversed(val))

	fields = ('chroms', 'genome_breaks', 'orients')
	new_vals = {}
	flip(fields, new_vals)

	if not change:
	    if not None in new_vals.values():
		returnable = []
		for field in fields:
		    returnable.append(new_vals[field])
		return returnable
	else:
	    if annotation:
		annot_fields = ('transcripts', 'transcript_breaks', 'exons', 'exon_bounds')
		flip(annot_fields, new_vals)
	    if not None in new_vals.values():
		for field, val in new_vals.iteritems():
		    setattr(self, field, val)
	            
    def key(self):
	"""Constructs a unique key for grouping adjacencies
	
	Args:
	    transcriptome: (boolean) whether adjacency is genomic or transcriptomic
	Returns:
	    A string that is used for grouping adjacencies
	"""
	if self.event is not None:
	    info = [self.event]
	else:
	    info = [self.rearrangement]

	chroms = self.chroms
	genome_breaks = self.genome_breaks
	orients = self.orients
	if not self.is_genome_breaks_sorted():
	    chroms, genome_breaks, orients = self.reverse_genome_breaks()

	for i in (0,1):
	    info.append(chroms[i])
	    info.append(genome_breaks[i])
	    if self.orients and\
	       ((type(orients) is tuple or type(orients) is list) and\
	        len(orients) == 2):
		info.append(orients[i])
	    else:
		info.append('na')
	return '-'.join(map(str, info))
    
    def update_transcript(self, transcript):
	self.genome_breaks = [transcript.txt_coord_to_genome_coord(coord) for coord in self.breaks]
	self.exons = [transcript.txt_coord_to_exon(coord) for coord in self.breaks]
    
    def set_probe(self, query_seq, len_on_each_side=50):
	"""Sets probe sequence of adjacency

	Tries to put len_on_each_side bases on each side of breakpoint
	And tries to put minimum of 2xlen_on_each_side
	i.e. if one side is less than len_on_each_side, tries to compensate on the other side
	"""
	breaks = sorted(self.seq_breaks)

	up_seq = query_seq[:breaks[0]]
	down_seq = query_seq[breaks[1] - 1:]
	mid_seq = query_seq[breaks[0]:breaks[1] - 1]
	up_short = down_short = 0
	if len(up_seq) < len_on_each_side:
	    up_short = len_on_each_side - len(up_seq)
	if len(down_seq) < len_on_each_side:
	    down_short = len_on_each_side - len(down_seq)
	if up_short > 0 and down_short > 0:
	    self.probe = query_seq
	elif up_short > 0:
	    self.probe = up_seq + mid_seq + down_seq[:min(len_on_each_side + up_short, len(down_seq))]
	elif down_short > 0:
	    self.probe = up_seq[-1 * min(len(up_seq), len_on_each_side + down_short)::] + mid_seq + down_seq
	else:
	    self.probe = up_seq[-1 * min(len(up_seq), len_on_each_side)::] + mid_seq + down_seq[:min(len_on_each_side, len(down_seq))]
	
    def get_subseqs(self, query_seq, len_on_each_side=None, seq_breaks=None):
	if not seq_breaks:
	    if type(self.seq_breaks) is tuple:
		breaks = map(int, sorted(list(self.seq_breaks)))
	else:
	    breaks = map(int, sorted(list(seq_breaks)))
	subseqs = []
	if breaks:
	    if len_on_each_side is None:
		subseqs.append(query_seq[:breaks[0]])
		subseqs.append(query_seq[breaks[1] - 1:])
	    else:
		subseqs.append(query_seq[max(0, breaks[0] - len_on_each_side) : breaks[0]])
		subseqs.append(query_seq[breaks[1] - 1 : min(len(query_seq), breaks[1] + len_on_each_side - 1)])

	return subseqs
	
    def extract_subseqs(self, contig_fasta):
	subseqs = []
	aligns = self.aligns[0]
	for i in range(len(aligns)):
	    subseqs.append(contig_fasta.fetch(self.contigs[0], aligns[i].qstart - 1, aligns[i].qend))
	
	return subseqs
    		        
    @classmethod
    def merge(cls, all_adjs):
	"""Merge adjacencies that have the same breakpoint (and same event type) together
	Args:
	    adjs: (list) Adjacency
	    transcriptome: (boolean) whether adjacency is genomic or transcriptomic
	Returns:
	    List of adjs with subsets that represent the same adjacency merged
	"""
	def merge_attr(adjs, attr):
	    values = []
	    for adj in adjs:
		value = getattr(adj, attr)
		if type(value) is tuple or type(value) is list:
		    value = '%s-%s' % (value[0], value[1]) 
		values.append(value)
	    return ','.join(map(str, values))

	merged = defaultdict(list)
	for adj in all_adjs:
	    adj.update_support_span()
	    merged[adj.key()].append(adj)
	    
	merged_adjs = []
	for key, adjs in merged.iteritems():
	    if len(adjs) == 1:
		for attr in ('seq_id', 'seq_breaks', 'support_span', 'probe'):
		    setattr(adjs[0], attr, merge_attr([adjs[0]], attr))
		merged_adjs.append(adjs[0])
	    else:
		merged_adj = deepcopy(adjs[0])
		for attr in ('seq_id', 'seq_breaks', 'support_span', 'probe'):
		    setattr(merged_adj, attr, merge_attr(adjs, attr))
		merged_adjs.append(merged_adj)
		
	return merged_adjs