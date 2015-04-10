from pybedtools import BedTool
from sets import Set

class Transcript:
    def __init__(self, id, gene=None, strand=None, coding=False):
	self.id = id
	self.gene = gene
	self.strand = strand
	self.exons = []
	self.coding = coding
	self.cds_start = None
	self.cds_end = None
	
    def add_exon(self, exon):
	self.exons.append(exon)
	self.exons.sort(key=lambda e: int(e[0]))
	    
    def exon(self, num):
	assert type(num) is int, 'exon number %s not given in int' % num
	assert self.strand == '+' or self.strand == '-', 'transcript strand not valid: %s %s' % (self.id, self.strand)
	
	if num < 1 or num > len(self.exons):
	    print 'exon number out of range:%s (1-%d)' % (num, len(self.exons))
	    return None
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
		#if feature.attrs.has_key('exon_number'):
		    #exon_num = int(feature.attrs['exon_number'])
		transcript_id = feature.attrs['transcript_id']
		gene = None
		if feature.attrs.has_key('gene_name'):
		    gene = feature.attrs['gene_name']
		elif feature.attrs.has_key('gene_id'):
		    gene = feature.attrs['gene_id']
		strand = feature.strand
		
		if feature.attrs.has_key('gene_biotype') and feature.attrs['gene_biotype'] == 'protein_coding':
		    coding = True
		else:
		    coding = False
				
		try:
		    transcript = transcripts[transcript_id]
		except:
		    transcript = Transcript(transcript_id, gene=gene, strand=strand, coding=coding)
		    transcripts[transcript_id] = transcript
		    
		transcript.add_exon(exon)
		
	    elif feature[2] == 'CDS':
		transcript_id = feature.attrs['transcript_id']
		strand = feature.strand
		cds = (int(feature.start) + 1, int(feature.stop))
		
		try:
		    transcript = transcripts[transcript_id]
		except:
		    transcript = Transcript(transcript_id, gene=gene, strand=strand, coding=coding)
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
