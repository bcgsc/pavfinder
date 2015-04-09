from sets import Set
import re
import gzip
from pavfinder.shared.adjacency import Adjacency
#from pavfinder.splice.fusion_finder import screen_realigns

class Event(Adjacency):
    # headers of tab-delimited output
    headers = ['ID',
               'event',
               'chrom1',
               'pos1',
               'orient1',
               'chrom2',
               'pos2',
               'orient2',
               'size',
               'contigs',
               'contig_breaks',
               'contig_support_span',
               'homol_seq',
               'homol_coords',
               'homol_len',
               'novel_sequence',
               'gene1',
               'transcript1',
               'exon1',
               'exon_bound1',
               'gene2',
               'transcript2',
               'exon2',
               'exon_bound2',
               'sense_fusion',
               "5'gene",
               "3'gene",
               'support_reads',
               'jn_depth',
               'ref5_jn_coord',
               'ref5_jn_depth',
               'ref3_jn_coord',
               'ref3_jn_depth',
               ]

    event_types = ('novel_donor',
                   'novel_acceptor',
                   'novel_exon',
                   'skipped_exon',
                   'novel_intron',
                   'retained_intron')
    
    def __init__(self, *args, **kwargs):
	Adjacency.__init__(self, *args, **kwargs)
	self.read_depth = 0
	self.ref5_coord = None
	self.ref5_depth = None
	self.ref3_coord = None
	self.ref3_depth = None
    
    @classmethod
    def output(cls, events, outdir, sort_by_event_type=False):
	"""Output events

	Args:
	    events: (list) Adjacency
	    outdir: (str) absolute path of output directory
	Returns:
	    events will be output in file outdir/events.tsv
	"""
	def get_smaller_pos(event):
	    """Returns 'smaller' coordinate of given event"""
	    if cls.compare_pos((event.chroms[0], event.breaks[0]), (event.chroms[1], event.breaks[1])) > 0:
		return (event.chroms[1], event.breaks[1])
	    else:
		return (event.chroms[0], event.breaks[0])
	    
	event_handlings = {
	    'fusion': 'from_fusion',
	    'ITD': 'from_single_locus',
	    'PTD': 'from_single_locus',
	    'ins': 'from_single_locus',
	    'del': 'from_single_locus',
	    'dup': 'from_single_locus',
	    'inv': 'from_single_locus',
	    'skipped_exon': 'from_single_locus',
	    'novel_exon': 'from_single_locus',
	    'novel_donor': 'from_single_locus',
	    'novel_acceptor': 'from_single_locus',
	    'novel_intron': 'from_single_locus',
	    'retained_intron': 'from_single_locus',
	}

	out_file = '%s/events.tsv' % outdir
	out = open(out_file, 'w')
	out.write('%s\n' % '\t'.join(cls.headers))
	
	if sort_by_event_type:
	    events_sorted = []
	    event_types = ['fusion',
	                   'ITD',
	                   'PTD',
	                   'ins',
	                   'del',
	                   'dup',
	                   'inv',
	                   'skipped_exon',
	                   'novel_exon',
	                   'novel_donor',
	                   'novel_acceptor',
	                   'novel_intron',
	                   'retained_intron',
	                   ]
	    for event_type in event_types:
		events_sorted.extend(sorted([e for e in events if e.rna_event == event_type], 
		                            cmp=lambda x,y: cls.compare_pos(get_smaller_pos(x), get_smaller_pos(y))))
	
	else:
	    events_sorted = sorted(events, cmp=lambda x,y: cls.compare_pos(get_smaller_pos(x), get_smaller_pos(y)))
	
	for event in events_sorted:
	    if event.rna_event:
		out_line = getattr(event, event_handlings[event.rna_event])()
		if out_line:
		    out.write('%s\n' % out_line)
	out.close()
	
    @staticmethod
    def output_reads(events, support_reads, outfile):
	"""Outputs support spanning reads for events"""
	recs = ''
	for event in events:
	    key = event.key(transcriptome=True)
	    if support_reads.has_key(key):
		for read_name, seq in support_reads[key]:
		    recs += '>%s-%s\n%s\n' % (key, read_name, seq)

	gzipped_outfile = outfile + '.gz'
	out = gzip.open(gzipped_outfile, 'wb')
	out.write(recs)
	out.close()
		    
    def from_fusion(self):
	"""Generates output line for a fusion/PTD event"""
	data = [self.id, self.rna_event]
	
	# sort breakpoints for output
	paired_values = []
	for values in zip(self.chroms, self.breaks, self.orients, self.genes, self.transcripts, self.exons, self.exon_bound):
	    paired_values.append(values)
	if Event.compare_pos((self.chroms[0], self.breaks[0]), (self.chroms[1], self.breaks[1])) > 0:
	    paired_values.reverse()
	for values in paired_values:
	    data.extend(values[:3])
	    
	# size not applicable to fusion
	data.append('-')
	
	# contigs and contig breaks and contig support span
	data.append(','.join(self.contigs))
	data.append(Event.to_string(self.contig_breaks))
	data.append(Event.to_string(self.contig_support_span))
		    
	# homol_seq and coords
	if self.homol_seq:
	    data.append(self.homol_seq[0])
	else:
	    data.append('-')
	if self.homol_coords:
	    homol_coords = []
	    for coords in self.homol_coords:
		homol_coords.append('-'.join(map(str, coords)))
	    data.append(';'.join(homol_coords))
	else:
	    data.append('-')
	if self.homol_seq:
	    data.append(len(self.homol_seq[0]))
	else:
	    data.append('-')
		    
	# novel_seq
	if hasattr(self, 'novel_seq') and self.novel_seq is not None:
	    data.append(self.novel_seq)
	else:
	    data.append('-')
	
	# gene, transcripts, exons, exon_bounds
	for values in paired_values:
	    data.extend(values[3:])
	
	# sense fusion, 5'gene, 3'gene
	data.append(self.is_sense)
	data.append(self.gene5)
	data.append(self.gene3)
	
	# support
	if not self.support or not self.support['spanning']:
	    data.append('-')
	else:
	    data.append(self.support['spanning'])

	# ref and alt depths
	for attr in ('read_depth', 'ref5_coord', 'ref5_depth', 'ref3_coord', 'ref3_depth'):
	    value = None
	    if hasattr(self, attr):
		value = getattr(self, attr)

	    if value is not None:
		data.append(getattr(self, attr))
	    else:
		data.append('-')
	
	return '\t'.join(map(str, data))

    def from_single_locus(self):
	"""Generates output line for an event from a single alignment
	
	Args:
	    event: (Adjacency) indel, ITD, splicing event (coming from single alignment)
	Returns:
	    Tab-delimited line
	"""
	data = [self.id, self.rna_event]
	
	chroms = (self.chroms[0], self.chroms[0])
	orients = ('L', 'R')
	for values in zip(chroms, self.breaks, orients):
	    data.extend(values)
	    
	# size
	if hasattr(self, 'size') and self.size is not None:
	    data.append(self.size)
	else:
	    data.append('-')
	    
	# contigs and contig breaks and contig support span
	data.append(','.join(self.contigs))
	data.append(Event.to_string(self.contig_breaks))
	data.append(Event.to_string(self.contig_support_span))
	
	# homol_seq and coords
	data.append('-')
	data.append('-')
	data.append('-')
	
	# novel seq
	if hasattr(self, 'novel_seq') and self.novel_seq is not None:
	    data.append(self.novel_seq)
	else:
	    data.append('-')

	# gene, transcripts, exons, exon_bounds
	genes = (self.genes[0], self.genes[0])
	transcripts = (self.transcripts[0], self.transcripts[0])
	if self.exons:
	    if len(self.exons) == 2:
		exons = self.exons
	    else:
		exons = (self.exons[0], self.exons[0])
	# novel exons
	else:
	    exons = ('-', '-')
	    
	# exon_bound
	exon_bound = ('-', '-')
	
	for values in zip(genes, transcripts, exons, exon_bound):
	    data.extend(values)
	    	    
	# sense fusion, 5'gene, 3'gene
	data.append('-')
	data.append('-')
	data.append('-')

	#support
	if not self.support or not self.support['spanning']:
	    data.append('-')
	else:
	    data.append(self.support['spanning'])

	# ref and alt depths
	for attr in ('read_depth', 'ref5_coord', 'ref5_depth', 'ref3_coord', 'ref3_depth'):
	    value = None
	    if hasattr(self, attr):
		value = getattr(self, attr)

	    if value is not None:
		data.append(getattr(self, attr))
	    else:
		data.append('-')

	return '\t'.join(map(str, data))
	
    @staticmethod
    def to_string(value):
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
	    
    @staticmethod
    def filter_by_support(events, min_support):
	"""Filters out events that don't have minimum spanning read support
	
	Args:
	    events: (list) Event
	    min_support: (int) minimum spanning read support
	"""
	out_indices = []
	for i in reversed(range(len(events))):
	    if not events[i].support['spanning'] or events[i].support['spanning'] < min_support:
		out_indices.append(i)
		
	for i in out_indices:
	    del events[i]
	    
    @staticmethod
    def compare_pos(pos1, pos2):
	"""Compares 2 genomic positions
	
	Args:
	    pos1: (tuple) chromosome1, coordinate1
	    pos2: (tuple) chromosome2, coordinate2
	Returns:
	    1 : pos2 > pos1
	    -1: pos1 < pos2
	    0 : pos1 == pos2
	"""
	chr1, coord1 = pos1
	chr2, coord2 = pos2
	
	if chr1[:3].lower() == 'chr':
	    chr1 = chr1[3:]
	if chr2[:3].lower() == 'chr':
	    chr2 = chr2[3:]
	
	if re.match('^\d+$', chr1) and not re.match('^\d+$', chr2):
	    return -1
	elif not re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
	    return 1
	else:
	    if re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
		chr1 = int(chr1)
		chr2 = int(chr2)
		
	    if chr1 < chr2:
		return -1
	    elif chr1 > chr2:
		return 1
	    else:
		if int(coord1) < int(coord2):
		    return -1
		elif int(coord1) > int(coord2):
		    return 1
		else:
		    return 0

    def is_splicing_event(self):
	return self.rna_event in self.event_types