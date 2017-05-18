from pybedtools import BedTool
from sets import Set
from alignment import reverse_complement
from translate import get_orfs
from math import ceil

class Transcript:
    def __init__(self, id, gene=None, strand=None, coding=False, chrom=None):
	self.id = id
	self.gene = gene
	self.strand = strand
	self.exons = []
	self.coding = coding
	self.cds_start = None
	self.cds_end = None
	self.chrom = chrom
	
    def add_exon(self, exon):
	self.exons.append(exon)
	self.exons.sort(key=lambda e: int(e[0]))
	    
    def exon(self, num, transcript_coord=False):
	assert type(num) is int, 'exon number %s not given in int' % num
	assert self.strand == '+' or self.strand == '-', 'transcript strand not valid: %s %s' % (self.id, self.strand)
	
	if num < 1 or num > len(self.exons):
	    print 'exon number out of range:%s (1-%d)' % (num, len(self.exons))
	    return None

	if not transcript_coord:
	    if self.strand == '+':
		return self.exons[num - 1]
	    else:
		return self.exons[len(self.exons) - num]
	else:
	    return self.exons_in_transcript_coords()[num - 1]

    def exons_in_transcript_coords(self):
	"""Returns exons in transcript coordinates"""
	if self.strand == '+':
	    blocks = self.exons
	else:
	    blocks = self.exons[::-1]

	exons = []
	start = 1
	for block in blocks:
	    exons.append([start, start + block[1] - block[0]])
	    start = exons[-1][1] + 1
	return exons
	
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
	
    def translate(self, fasta):
	if self.is_coding():
	    cdna = self.get_sequence(fasta, cds_only=True)
	    orfs = get_orfs(cdna, frames=[0], only_strand='+')
	    if orfs:
		return orfs[0][-1]

	return None

    def txt_coord_to_aa_coord(self, txt_coord):
	return int(ceil(float(txt_coord) / 3))

    def aa_coord_to_genome_coord(self, aa_coord, cdna_seq=None):
	if self.is_coding() and self.cds_start is not None and self.cds_end is not None:
	    cds_coords = (aa_coord - 1) * 3 + 1, (aa_coord - 1) * 3 + 3
	    utr5 = 0
	    if self.strand == '+':
		for i in range(len(self.exons)):
		    exon = self.exons[i]
		    if exon[1] < self.cds_start:
			utr5 += self.get_exon_size(i)
		    elif self.cds_start >= exon[0] and self.cds_start <= exon[1]:
			utr5 += self.cds_start - exon[0]
			break
	    else:
		for i in range(len(self.exons))[::-1]:
		    exon = self.exons[i]
		    if self.cds_start < exon[0]:
			utr5 += self.get_exon_size(i)
		    elif self.cds_start >= exon[0] and self.cds_start <= exon[1]:
			utr5 += exon[1] - self.cds_start
			break
	    txt_coords = utr5 + cds_coords[0], utr5 + cds_coords[1]
	    return self.txt_coord_to_genome_coord(txt_coords[0]),\
	           self.txt_coord_to_genome_coord(txt_coords[1])
	return None

    def aa_coord_to_txt_coord(self, aa_coord, cds=True):
	if self.is_coding():
	    if cds:
		return (aa_coord - 1) * 3 + 1, (aa_coord - 1) * 3 + 3
	return None

    def get_sequence(self, fasta, cds_only=False, genomic=False):
	"""Extract transcript sequence using CDS coordinate info

	if "cds_only", sequence will begin with start codon ATG and
	end with nucleotide before stop codon

	Args:
	    fasta - Pysam FastaFile
	    cds_only(optional) - coding sequence only
	Returns:
	    sequence will be same strand of transcript
	    i.e. reverse_complemented if strand of transcript is negative
	    None if cds asked but transcript object doesn't have cds start and end coordinates
	"""
	seq = ''
	inside_cds = False

	exons = self.exons
	if self.strand == '-':
	    exons = reversed(exons)

	if cds_only and (self.cds_start is None or self.cds_end is None):
	    return None
	
	if genomic:
	    return fasta.fetch(self.chrom, self.exons[0][0] - 1, self.exons[-1][1])

	for exon in exons:
	    span = None
	    if not cds_only:
		span = exon[0] - 1, exon[1]
	    else:
		if not inside_cds and self.cds_start >= exon[0] and self.cds_start <= exon[1]:
		    if self.strand == '+':
			span = self.cds_start - 1, exon[1]
		    else:
			span = exon[0] - 1, self.cds_start
		    inside_cds = True
		elif inside_cds and self.cds_end >= exon[0] and self.cds_end <= exon[1]:
		    if self.strand == '+':
			span = exon[0] - 1, self.cds_end
		    else:
			span = self.cds_end - 1, exon[1]
		    inside_cds = True
		elif cds_only and inside_cds:
		    span = exon[0] - 1, exon[1]

	    if span is not None:
		if span is not None:
		    exon_seq = fasta.fetch(self.chrom, span[0], span[1])
		if self.strand == '-':
		    exon_seq = reverse_complement(exon_seq)
		seq += exon_seq.upper()

	return seq

    def get_exon_size(self, index):
	"""Returns exon size given exon index (not exon number)

	index = {0,1,...,len(exons)-1}, already from left to right coordinates
	"""
	if index >= 0 and index < len(self.exons):
	    return self.exons[index][1] - self.exons[index][0] + 1
	return None

    def coord_to_exon(self, coord):
	"""Given genomic coordinate, return exon number"""
	for i in range(len(self.exons)):
	    if coord >= self.exons[i][0] and coord <= self.exons[i][1]:
		return self.exon_num(i)

    def txt_coord_to_exon(self, coord):
	"""Given transcript coordinate, return exon number

	Motivation: For converting contig to transrcipt alignment back to
	            genomic coordinate reporting
	"""
	exon_start = 1
	indexes = range(len(self.exons))
	if self.strand == '-':
	    indexes.reverse()
	for i in range(len(indexes)):
	    if i > 0:
		exon_start += self.get_exon_size(indexes[i - 1])
	    if coord >= exon_start and coord <= exon_start + self.get_exon_size(indexes[i]) - 1:
		return self.exon_num(indexes[i])
	return None
    
    def txt_genomic_coord_to_genome_coord(self, coord):
	return self.exons[0][0] + coord - 1

    def txt_coord_to_genome_coord(self, coord):
	"""Given transcript coordinate, returns corresponding genomic coordinate"""
	exon_start = 1
	indexes = range(len(self.exons))
	if self.strand == '-':
	    indexes.reverse()
	for i in range(len(indexes)):
	    if i > 0:
		exon_start += self.get_exon_size(indexes[i - 1])
	    if coord >= exon_start and coord <= exon_start + self.get_exon_size(indexes[i]) - 1:
		if self.strand == '+':
		    return self.exons[indexes[i]][0] + (coord - exon_start)
		else:
		    return self.exons[indexes[i]][1] - (coord - exon_start)
	return None
    
    def genome_coord_to_txt_coord(self, coord, cds=False):
	"""Given genome coordinate, returns corresponding transcript coordinate"""
	def get_exon_size(exons, index):
	    return exons[index][1] - exons[index][0] + 1

	exons = self.exons
	if cds:
	    exons = self.get_cds()
	if exons:
	    for i in range(len(exons)):
		if coord >= exons[i][0] and coord <= exons[i][1]:
		    if self.strand == '+':
			offset = coord - exons[i][0] + 1
			tcoord = 0
			for j in range(i):
			    tcoord += get_exon_size(exons, j)
		    else:
			offset = exons[i][1] - coord + 1
			tcoord = 0
			for j in range(i + 1, len(exons)):
			    tcoord += get_exon_size(exons, j)

		    return tcoord + offset
	return None
    
    def get_cds(self):
	cds = []
	if self.is_coding() and self.cds_start is not None and self.cds_end is not None:
	    cds_start, cds_end = sorted([int(self.cds_start), int(self.cds_end)])

	    for exon in self.exons:
		if cds_start > exon[1]:
		    continue
		elif (cds_start >= exon[0] and cds_start <= exon[1]) or\
		     (cds_end >= exon[0] and cds_end <= exon[1]):
		    cds.append([max(cds_start, exon[0]), min(cds_end, exon[1])])
		elif cds_end < exon[0]:
		    break
		else:
		    cds.append(exon)
	return cds

    def at_exon_bound(self, genome_coord=None, txt_coord=None, exon_num=None, return_edge=False):
	# coord = genomic
	coord = None
	if genome_coord is not None:
	    coord = genome_coord
	elif txt_coord is not None:
	    coord = self.txt_coord_to_genome_coord(txt_coord)
	    
	if coord is not None:
	    if exon_num is not None:
		exon = self.exon(exon_num)
		if coord == exon[0] or coord == exon[1]:
		    if return_edge:
			if coord == exon[0]:
			    if self.strand == '+':
				return 5
			    else:
				return 3
			else:
			    if self.strand == '+':
				return 3
			    else:
				return 5
		    else:
			return True
	    else:
		for exon in self.exons:
		    if coord == exon[0] or coord == exon[1]:
			return True
	    return False
	return None
    
    def is_coding(self):
	if self.coding is None:
	    if self.cds_start is not None and self.cds_end is not None and\
	       self.cds_start != self.cds_end:
		self.coding = True
	return self.coding
	
    def utr(self, end=None):
	# cds_start and cds_end are in reference to transcript strand
	# i.e. cds_start > cds_end if transcript strand is '-'
	if self.is_coding and self.cds_start is not None and self.cds_end is not None:
	    if self.strand == '+':
		utr5 = self.exons[0][0], self.cds_start - 1
                utr3 = self.cds_end + 4, self.exons[-1][1]
            else:
                utr5 = self.cds_start + 1, self.exons[-1][1]
                utr3 = self.exons[0][0], self.cds_end - 4

	    if end == 5:
		return utr5
	    elif end == 3:
		return utr3
	    else:
		return utr5, utr3
	return None
	
    def within_utr(self, pos):
	utrs = self.utr()
	if utrs is not None:
	    utr5, utr3 = utrs
	    if pos >= utr5[0] and pos <= utr5[1]:
		return 5
	    elif pos >= utr3[0] and pos <= utr3[1]:
		return 3
	return False

    @staticmethod
    def extract_transcripts(annotation_file):
	"""Extracts all exon info into transcript objects
	
	Requires annotation file passed to object
	Uses PyBedTool for parsing
	
	Returns:
	    List of Transcripts with exon info, strand
	"""
	transcripts = {}
	for feature in BedTool(annotation_file):
	    if feature[2] == 'exon':
		exon = (int(feature.start) + 1, int(feature.stop))
		transcript_id = feature.attrs['transcript_id']
		gene = None
		if feature.attrs.has_key('gene_name'):
		    gene = feature.attrs['gene_name']
		elif feature.attrs.has_key('gene_id'):
		    gene = feature.attrs['gene_id']
		strand = feature.strand
		
		coding = None
		if feature.attrs.has_key('gene_biotype'):
		    if feature.attrs['gene_biotype'] == 'protein_coding':
			coding = True
				
		try:
		    transcript = transcripts[transcript_id]
		except:
		    transcript = Transcript(transcript_id, gene=gene, strand=strand, chrom=feature.chrom)
		    if coding is not None:
			transcript.coding = coding
		    transcripts[transcript_id] = transcript
		    
		transcript.add_exon(exon)
		
	    elif feature[2] == 'CDS':
		transcript_id = feature.attrs['transcript_id']
		gene = None
		if feature.attrs.has_key('gene_name'):
		    gene = feature.attrs['gene_name']
		elif feature.attrs.has_key('gene_id'):
		    gene = feature.attrs['gene_id']
		strand = feature.strand
		cds = (int(feature.start) + 1, int(feature.stop))
		coding = True
		
		try:
		    transcript = transcripts[transcript_id]
		    transcript.coding = True
		except:
		    transcript = Transcript(transcript_id, gene=gene, strand=strand, coding=coding, chrom=feature.chrom)
		    transcripts[transcript_id] = transcript
		    
		
		if strand == '+':
		    if transcript.cds_start is None or cds[0] < transcript.cds_start:
			transcript.cds_start = cds[0]
		    if transcript.cds_end is None or cds[1] > transcript.cds_end:
			transcript.cds_end = cds[1]
		else:
		    if transcript.cds_end is None or cds[0] < transcript.cds_end:
			transcript.cds_end = cds[0]
		    if transcript.cds_start is None or cds[1] > transcript.cds_start:
			transcript.cds_start = cds[1]
			
	for transcript in transcripts.values():
	    if not transcript.coding and\
	       transcript.cds_start is not None and\
	       transcript.cds_end is not None and\
	       transcript.cds_start != transcript.cds_end:
		transcript.coding = True
		
	return transcripts

    @staticmethod
    def extract_features(annotation_file):
	"""Extracts all exons and exon-junctions from given gtf file

	Arguments:
	    annotation_file: gtf file
	Returns:
	    {'exon': {('chr', start, end), ('chr', start, end)}
	     'junction': {('chr', start, end), ('chr', start, end)}
	"""
	features = {'exon': Set(), 'junction': Set()}

	transcripts = {}
	for exon in BedTool(annotation_file).filter(lambda f: f[2] == 'exon'):
	    coords = int(exon.start) + 1, int(exon.stop)
	    pos = exon.chrom, coords[0], coords[1]
	    features['exon'].add(pos)
	    try:
		exons = transcripts[exon.attrs['transcript_id']]
	    except:
		transcripts[exon.attrs['transcript_id']] = []
		exons = transcripts[exon.attrs['transcript_id']]
	    exons.append(pos)

	for exons in transcripts.values():
	    for i in range(len(exons) - 1):
		jn = (exons[i][0], exons[i][2], exons[i + 1][1])
		features['junction'].add(jn)

	return features
