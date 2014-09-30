from optparse import OptionParser
from itertools import groupby, chain
from pybedtools import create_interval_from_list, set_tempdir, BedTool
import pysam
import sys
import re
from shared import gmap
from shared.annotate import overlap
from shared.alignment import reverse_complement
from sets import Set
from intspan import intspan
from SV.split_align import call_event
from SV.variant import Adjacency

class Fusion:
    def __init__(self, gene1, gene2):
	self.gene1 = gene1
	self.gene2 = gene2
	self.contigs = []
	self.breaks = []
	self.txt1 = None
	self.txt2 = None
	self.exon1 = None
	self.exon2 = None
	self.junctions = None
	
    def as_tab(self):
	data = []
	data.append(self.gene1)
	data.append(self.gene2)
	data.append(','.join(self.contigs))
	
	return '\t'.join(data)

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
    
class Event:
    headers = ['event_type',
               'chrom1',
               'pos1',
               'gene1',
               'transcripts1',
               'exons1',
               'chrom2',
               'pos2',
               'gene2',
               'transcripts2',
               'exons2',
               'contigs',
               'contig_blocks',
               'contig_breaks',
               ]
    
    @classmethod
    def output(cls, events, outdir):
	# chimeras
	# svs
	#svs = [event for event in events if event.event_type == 'ins' or event.event_type == 'del']
	#cls.output_sv(svs, '%s/sv.tsv' % outdir)
	
	# splicing
	splicing_types = ['novel_exon', 'skipped_exon']
	splicing = [event for event in events if event.rna_event in splicing_types]
	cls.output_splicing(splicing, '%s/splicing.tsv' % outdir)
	
	# fusions
	fusions = [event for event in events if event.rna_event == 'fusion']
	cls.output_fusion(fusions, '%s/fusions.tsv' % outdir)
	
	# indels
	print len(events)
	events[0].id = '123'
	print 'abcd', events, events[0], events[0].rearrangement, events[0].breaks, events[0].chroms, events[0].genes, events[0].transcripts
	indels = [event for event in events if event.rearrangement in ('ins', 'del')]
	cls.output_indels(indels, '%s/indels.tsv' % outdir)
	
    @classmethod
    def output_sv(cls, events, outfile):
	out = open(outfile, 'w')
	out.write('%s\n' % '\t'.join(cls.headers))
	for event in events:
	    out.write('%s\n' % event.as_tab())
	out.close()
	
    @classmethod
    def output_indels(cls, events, outfile):
	def as_indel(event):
	    data = []
	    data.append(event.rna_event)
	    data.append(event.chroms[0])
	    data.append(event.breaks[0])
	    data.append(event.breaks[1])
	    data.append(event.size)
	    data.append(event.genes[0])
	    data.append(event.transcripts[0])
	    if event.exons:
		data.append(','.join(map(str, event.exons)))
	    else:
		data.append('-')
	    data.append(','.join(event.contigs))
	    data.append(cls.to_string(event.contig_breaks))
	    if event.novel_seq:
		data.append(event.novel_seq)
		
	    return '\t'.join(map(str, data))
	
	headers = ['event',
	           'chrom',
	           'pos1',
	           'pos2',
	           'size',
	           'gene',
	           'transcript',
	           'exon',
	           'contigs',
	           'contig_breaks',
	           'novel_sequence',
	           ]
	out = open(outfile, 'w')
	out.write('%s\n' % '\t'.join(headers))
	for event in events:
	    out.write('%s\n' % as_indel(event))
	out.close()
	
	
    @classmethod
    def output_splicing(cls, events, outfile):
	def as_splice_variant(event):
	    data = []
	    data.append(event.rna_event)
	    data.append(event.chroms[0])
	    data.append(','.join(map(str, event.breaks)))
	    data.append(event.genes[0])
	    data.append(event.transcripts[0])
	    if event.exons:
		data.append(','.join(map(str, event.exons)))
	    else:
		data.append('-')
	    data.append(','.join(event.contigs))
	    data.append(cls.to_string(event.contig_breaks))
	    
	    return '\t'.join(map(str, data))
	    
	headers = ['event',
	           'chrom',
	           'coords',
	           'gene',
	           'transcript',
	           'exon',
	           'contigs',
	           'contig_breaks',
	           ]
	out = open(outfile, 'w')
	out.write('%s\n' % '\t'.join(headers))
	for event in events:
	    print 'lll', event.rna_event, event.breaks, event.exons, event.contigs, event.contig_breaks, event.genes, event.transcripts
	    out.write('%s\n' % as_splice_variant(event))
	out.close()
	
    @classmethod
    def as_fusion(cls, event):
	data = []
	for values in zip(event.chroms, event.breaks, event.orients, event.genes, event.transcripts, event.exons):
	    print 'abc', values
	    data.extend(values)
	
	data.append(','.join(event.contigs))
	data.append(cls.to_string(event.contig_breaks))
	
	return '\t'.join(map(str, data))
	    
    @classmethod
    def output_fusion(cls, events, outfile):
	headers = ['chrom1',
	           'break1',
	           'gene1',
	           'transcript1',
	           'exon1',
	           'chrom2',
	           'break2',
	           'gene2',
	           'transcript2',
	           'exon2',
	           'contigs',
	           'contig_blocks',
	           'contig_breaks',
	           ]
	out = open(outfile, 'w')
	out.write('%s\n' % '\t'.join(headers))
	for event in events:
	    out.write('%s\n' % cls.as_fusion(event))
	out.close()
	
    @classmethod
    def to_string(cls, value):
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
    def is_ITD(cls, adj, align, contigs_fasta, shift_size=20):
	if adj.novel_seq:
	    contig_seq = contigs_fasta.fetch(adj.contigs[0])
	    print 'itd', contig_seq, align.strand, adj.novel_seq, len(adj.novel_seq)
	    
	    novel_seq = adj.novel_seq if align.strand == '+' else reverse_complement(adj.novel_seq)
	    matches = re.findall(novel_seq, contig_seq)
	    if len(matches) > 1:
		return True
	    else:
		size_novel_seq = len(adj.novel_seq)
		start = min(adj.contig_breaks) 
		without_ins_seq_original = contig_seq[:start] + contig_seq[start + size_novel_seq:]
		print align.query, 'itd xx', start, start + size_novel_seq, len(contig_seq[:start]), len(contig_seq[start + size_novel_seq:]), size_novel_seq
		for i in range(1, shift_size + 1):
		    without_ins_seq = contig_seq[:start - i] + contig_seq[start - i + size_novel_seq:]
		    if without_ins_seq_original==without_ins_seq:
			novel_seq = contig_seq[start - i: start - i + size_novel_seq]
			if len(re.findall(novel_seq, contig_seq)) > 1:
			    return True
		    without_ins_seq = contig_seq[:start + i] + contig_seq[start + i + size_novel_seq:]
		    if without_ins_seq_original==without_ins_seq:
			novel_seq = contig_seq[start + i : start + i + size_novel_seq]
			if len(re.findall(novel_seq, contig_seq)) > 1:
			    return True		
	    
	return False
	    
    
#class Event2:
    #headers = ['event_type',
               #'chrom1',
               #'pos1',
               #'gene1',
               #'transcripts1',
               #'exons1',
               #'chrom2',
               #'pos2',
               #'gene2',
               #'transcripts2',
               #'exons2',
               #'contigs',
               #'contig_blocks',
               #'contig_breaks',
               #]
    
    #def __init__(Adjacency, event_type, 
                 #chrom1=None, pos1=None, break1=None, gene1=None, transcripts1=[], exons1=[], 
                 #chrom2=None, pos2=None, break2=None, gene2=None, transcripts2=[], exons2=[], 
                 #contigs=[], contig_blocks=[], contig_breaks=[]):
	#self.event_type = event_type
	#self.chrom1 = chrom1
	#self.pos1 = pos1
	#self.gene1 = gene1
	#self.transcripts1 = transcripts1
	#self.exons1 = exons1
	#self.chrom2 = chrom2
	#self.pos2 = pos2
	#self.gene2 = gene2
	#self.transcripts2 = transcripts2
	#self.exons2 = exons2
	#self.contigs = contigs
	#self.contig_blocks = contig_blocks
	#self.contig_breaks = contig_breaks
	
    #def as_tab(self):
	#def to_string(value):
	    #if value is None:
		#return '-'
	    #elif type(value) is list or type(value) is tuple:
		#items = []
		#for item in value:
		    #if item is None:
			#items.append('-')
		    #elif type(item) is list or type(item) is tuple:
			#if item:
			    #items.append(','.join(map(str, item)))
		    #else:
			#items.append(str(item))
			
		#if items:
		    #return ';'.join(items)
		#else:
		    #return '-'
	    #else:
		#return str(value)
	
	
	#cols = []
	#for attr in self.headers:
	    #value = getattr(self, attr)
	    #cols.append(to_string(value))
		
	#return '\t'.join(cols)
    
    #@classmethod
    #def output(cls, events, outdir):
	## chimeras
	## svs
	#svs = [event for event in events if event.event_type == 'ins' or event.event_type == 'del']
	#cls.output_sv(svs, '%s/sv.tsv' % outdir)
	
	## splicing
	#splicing_types = ['novel_exon', 'skipped_exon']
	#splicing = [event for event in events if event.event_type in splicing_types]
	#cls.output_splicing(splicing, '%s/splicing.tsv' % outdir)
	
    #@classmethod
    #def output_sv(cls, events, outfile):
	#out = open(outfile, 'w')
	#out.write('%s\n' % '\t'.join(cls.headers))
	#for event in events:
	    #out.write('%s\n' % event.as_tab())
	#out.close()
	
    #@classmethod
    #def output_splicing(cls, events, outfile):
	#out = open(outfile, 'w')
	#out.write('%s\n' % '\t'.join(cls.headers))
	#for event in events:
	    #out.write('%s\n' % event.as_tab())
	#out.close()
	
class Mapping:
    """Mapping per alignment"""
    def __init__(self, contig, align_blocks, transcripts=[]):
	self.contig = contig
	self.transcripts = transcripts
	self.genes = list(Set([txt.gene for txt in self.transcripts]))
	self.align_blocks = align_blocks
	self.coverages = []
	
    def overlap(self):
	align_span = self.create_span(self.align_blocks)
	
	for transcript in self.transcripts:
	    exon_span = self.create_span(transcript.exons)
	    olap = exon_span.intersection(align_span)
	    #print 'olap', olap, len(olap), len(exon_span), float(len(olap)) / float(len(exon_span))
	    self.coverages.append(float(len(olap)) / float(len(exon_span)))
	    
    @classmethod
    def create_span(cls, blocks):
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
	scores = {}
	metrics = {}
	for transcript, matches in mappings:
	    metric = {'score': 0,
	              'from_edge': 0,
	              'txt_size': 0
	              }
	    #print 'choice', transcript.id, transcript.gene, transcript.strand, matches, align.tstart, align.tend, matches[0][0][0], matches[-1][-1][0]
	    
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
		    score += 1
		if matches[i][-1][1][1] == '=':
		    score += 1		
			
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
	"""Group by gene"""
	gene_mappings = []
	for gene, group in groupby(all_mappings, lambda m: m.genes[0]):
	    mappings = list(group)
	    contigs = ','.join([mapping.contig for mapping in mappings])
	    transcripts = [mapping.transcripts for mapping in mappings]
	    #print transcripts
	    #print 'tt', Set(chain(*transcripts))
	    align_blocks = [mapping.align_blocks for mapping in mappings]
	    #print align_blocks
	    #print list(chain(*align_blocks))
	    
	    align_blocks = None
	    
	    for mapping in mappings:
		try:
		    align_blocks = align_blocks.union(cls.create_span(mapping.align_blocks))
		except:
		    align_blocks = cls.create_span(mapping.align_blocks)
	    #print align_blocks.ranges()
		
	    
	    gene_mappings.append(Mapping(contigs,
	                                 align_blocks.ranges(),
	                                 list(Set(chain(*transcripts))),
	                                 )
	                         )
	    
	#for mapping in gene_mappings:
	    #mapping.overlap()
	    #print 'gene', mapping.as_tab()
	    
	[mapping.overlap() for mapping in gene_mappings]
	return gene_mappings
    
    @classmethod
    def output(cls, mappings, outfile):
	out = open(outfile, 'w')
	out.write('%s\n' % cls.header())
	for mapping in mappings:
	    out.write('%s\n' % mapping.as_tab())
	out.close()

class ExonMapper:
    def __init__(self, bam_file, aligner, contigs_fasta_file, annotation_file, ref_fasta_file, outdir):
        self.bam = pysam.Samfile(bam_file, 'rb')
	self.contigs_fasta = pysam.Fastafile(contigs_fasta_file)
        self.ref_fasta = pysam.Fastafile(ref_fasta_file)
	self.annot = pysam.Tabixfile(annotation_file, parser=pysam.asGTF())
        self.aligner = aligner
        self.outdir = outdir
        
        self.blocks_bed = '%s/blocks.bed' % outdir
        self.overlaps_bed = '%s/blocks_olap.bed' % outdir
        self.aligns = {}        

        self.annotations_file = annotation_file
	
	self.mappings = []
	self.events = []
	
	transcripts = self.extract_transcripts()
	self.map_contigs_to_transcripts(transcripts)
		
    def get_attrs(self, attrs_str):
	attrs = {}
	for field in attrs_str.split(';'):
	    cols = field.rstrip(' ').lstrip(' ').split(' ')
	    if len(cols) == 2:
		attrs[cols[0]] = cols[1].strip('"')
	    
	return attrs
	
    def map_contigs_to_transcripts(self, transcripts):
	aligns = []
	for contig, group in groupby(self.bam.fetch(until_eof=True), lambda x: x.qname):
            alns = list(group)
	    mappings = []
	    
	    aligns = self.extract_aligns(alns)
	    if aligns is None:
		print 'no valid alignment', contig
		continue
	    
	    chimera = True if len(aligns) > 1 else False
	    chimera_block_matches = []
	    	    
	    for align in aligns:
		if align is None:
		    print 'bad alignment', contig
		    continue
		
		if '_' in align.target:
		    print 'bad target', contig, align.target
		    continue
		
		print 'bbb', align.query, align.blocks, align.query_blocks
		
		# entire contig align within single exon or intron
		within_intron = []
		within_exon = []
		
		transcripts_mapped = Set()
		if re.search('[._]', align.target):
		    continue
		
		for gtf in self.annot.fetch(align.target, align.tstart, align.tend):			
		    if not chimera and align.tstart >= gtf.start and align.tend <= gtf.end:
			match = self.match_exon((align.tstart, align.tend), (gtf.start, gtf.end)) 
			
			if gtf.feature == 'intron':	
			    within_intron.append((gtf, match))
			    
			elif gtf.feature == 'exon':
			    within_exon.append((gtf, match))
			       
		    else:
			attrs = self.get_attrs(gtf.attributes)
			if gtf.feature == 'exon':
			    transcripts_mapped.add(gtf.transcript_id)
			    		
		if transcripts_mapped:
		    # key = transcript name, value = "matches"
		    all_block_matches = {}
		    # Transcript objects that are fully matched
		    full_matched_transcripts = []
		    for txt in transcripts_mapped:
			#print 'gene', contig, len(align.blocks), transcripts[txt].id, transcripts[txt].gene, transcripts[txt].strand, transcripts[txt].num_exons()			
			block_matches = self.map_exons(align.blocks, transcripts[txt].exons)
			all_block_matches[txt] = block_matches
			
			if not chimera and self.is_full_match(block_matches):
			    print 'full', contig
			    full_matched_transcripts.append(transcripts[txt])
			mappings.append((transcripts[txt], block_matches))
			    
		    if not full_matched_transcripts:			
			events = self.find_events(all_block_matches, align, transcripts)
			if not events:
			    print 'partial_no_events', contig
			else:
			    print 'events', align.query, events
			    self.events.extend(events)
			    			    
		    if chimera:
			chimera_block_matches.append(all_block_matches)
		    
		elif not chimera:
		    if within_exon:
			print 'exon', contig, align.target, align.tstart, align.tend#, within_exon[0].gene_id, within_exon[0].start, within_exon[0].end
		    
		    elif within_intron:
			#print 'aa', within_intron
			print 'intron', contig, align.target, align.tstart, align.tend#, within_intron[0].gene_id, within_intron[0].start, within_intron[0].end
		
		#if chimera:
		    #chimera_block_matches.append(all_block_matches)
		
	    if chimera and chimera_block_matches:
		if len(chimera_block_matches) == len(aligns):
		    print 'chimera', chimera_block_matches
		    fusion = self.find_chimera(chimera_block_matches, transcripts, aligns)
		    if fusion:
			self.events.append(fusion)
		
	   
	    # pick best matches
	    if len(mappings) > 1:
		best_mapping = Mapping.pick_best(mappings, align, debug=True)
		self.mappings.append(best_mapping)
		
	
    def map_exons(self, blocks, exons):
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
		    
    #def extract_annot(self):
	#self.annot = {}
	#for feature in BedTool(self.annotations_file):
	    #if feature[2] == 'exon':		
		#print feature
		#try:
		    #self.annot[feature.attrs['transcript_id']]['exons'][int(feature.attrs['exon_number'])] = (int(feature.start) + 1, int(feature.stop))
		#except:
		    #self.annot[feature.attrs['transcript_id']] = {}
		    #self.annot[feature.attrs['transcript_id']]['gene'] = feature.attrs['gene_name']
		    #self.annot[feature.attrs['transcript_id']]['strand'] = feature.strand
		    #self.annot[feature.attrs['transcript_id']]['exons'] = {}
		    #self.annot[feature.attrs['transcript_id']]['exons'][int(feature.attrs['exon_number'])] = (int(feature.start) + 1, int(feature.stop))
		    
    def extract_transcripts(self):
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
		
	#for transcript in transcripts.values():
	    #print transcript.id, transcript.num_exons(), transcript.strand
	    
	return transcripts
    
    def identify_fusion(self, matches1, matches2, transcripts):
	"""Given 2 block matches pick the 2 transcripts"""
	print 'iii', matches1, matches2
	scores = {}
	for transcript in matches1:
	    score = 0
	    if matches1[transcript] is not None:
		if matches1[transcript][0][1][1] == '=':
		    score += 100
		    
		# align block within exon gets more points
		if matches1[transcript][0][1][0] == '=':
		    score += 15
		elif matches1[transcript][0][1][0] == '>':
		    score += 10
		elif matches1[transcript][0][1][0] == '<':
		    score += 5
	    scores[transcript] = score
	    
	best_score = max(scores.values())
	#print best_score, [t for t in matches1.keys() if scores[t] == best_score]
	best_txt1 = sorted([t for t in matches1.keys() if scores[t] == best_score], 
	                   key=lambda t: transcripts[t].length(), reverse=True)[0]
	#print best_txt1
	
	
	scores = {}
	for transcript in matches2:
	    score = 0
	    if matches2[transcript] is not None:
		if matches2[transcript][0][1][1] == '=':
		    score += 100
		    
		if matches2[transcript][0][1][0] == '=':
		    score += 15
		elif matches2[transcript][0][1][0] == '>':
		    score += 10
		elif matches2[transcript][0][1][0] == '<':
		    score += 5
	    scores[transcript] = score
	    
	best_score = max(scores.values())
	#print best_score, [t for t in matches2.keys() if scores[t] == best_score]
	best_txt2 = sorted([t for t in matches2.keys() if scores[t] == best_score], 
	                   key=lambda t: transcripts[t].length(), reverse=True)[0]
	#print best_txt2
	print 'fusion', best_txt1, matches1[best_txt1], best_txt2, matches2[best_txt2]
	return (best_txt1, matches1[best_txt1]), (best_txt2, matches2[best_txt2])
    
    def identify_fusion_unknown_break(self, matches, transcripts):
	all_matches = Set()
	for txt in matches:
	    if matches[txt] is None:
		continue
	    
	    for match in matches[txt]:
		all_matches.add((txt, match))
		
	left_bound_matches = [match for match in all_matches if match[1][1][0] == '=']
	right_bound_matches = [match for match in all_matches if match[1][1][1] == '=']

	print 'left', left_bound_matches
	print 'right', right_bound_matches
	
	if left_bound_matches and right_bound_matches:
	    left_bound_matches.sort(key=lambda m: transcripts[m[0]].length(), reverse=True)
	    right_bound_matches.sort(key=lambda m: transcripts[m[0]].length(), reverse=True)
	    
	    print 'left', left_bound_matches
	    print 'right', right_bound_matches
	    
	    print 'fusion_unknown_break', left_bound_matches[0], right_bound_matches[0]
	    
    def find_chimera(self, chimera_block_matches, transcripts, aligns):	
	assert len(chimera_block_matches) == len(aligns), 'number of block matches(%d) != number of aligns(%d)' % \
	       (len(chimera_block_matches), len(aligns))
	for (matches1_by_txt, matches2_by_txt) in zip(chimera_block_matches, chimera_block_matches[1:]):
	    print 'a', matches1_by_txt
	    print 'b', matches2_by_txt
	    
	    genes1 = Set([transcripts[txt].gene for txt in matches1_by_txt.keys()])
	    genes2 = Set([transcripts[txt].gene for txt in matches2_by_txt.keys()])
	    	    
	    print genes1, genes2
	    # gene fusion
	    print 'g1', genes1
	    print 'g2', genes2
	    
	    junc_matches1 = {}
	    num_blocks = len(aligns[0].blocks)
	    for transcript in chimera_block_matches[0].keys():
		#print 'ggg', transcript, chimera_block_matches[0][transcript]
		junc_matches1[transcript] = chimera_block_matches[0][transcript][num_blocks - 1]
	
	    junc_matches2 = {}
	    for transcript in chimera_block_matches[1].keys():
		junc_matches2[transcript] = chimera_block_matches[1][transcript][0]
    
	    print 'j1', junc_matches1
	    print 'j2', junc_matches2
    
	    #print 'contig', aligns[0].query
	    junc1, junc2 = self.identify_fusion(junc_matches1, junc_matches2, transcripts)
	    print 'fusion-chimera', aligns[0].query, transcripts[junc1[0]].gene, transcripts[junc2[0]].gene
	    
	    if junc1 and junc2:
		fusion = call_event(aligns[0], aligns[1])
		fusion.rna_event = 'fusion'
		#print 'ggg', event.event_type, event.chroms, event.breaks, event.rearrangement, event.contig_breaks, event.orients
		#print 'kkk0', adj.chroms, adj.breaks, adj.rearrangement, adj.contig_breaks, adj.orients
		print 'kkk', junc1, junc2, transcripts[junc1[0]], transcripts[junc2[0]], fusion.breaks, fusion.contig_breaks
		print 'kkk1', transcripts[junc1[0]].gene, junc1[0], junc1[1][0][0], transcripts[junc1[0]].strand, aligns[0].target, aligns[0].qstart, aligns[0].qend, aligns[0].strand
		print 'kkk2', transcripts[junc2[0]].gene, junc2[0], junc2[1][0][0], transcripts[junc2[0]].strand, aligns[1].target, aligns[1].qstart, aligns[1].qend, aligns[1].strand

		fusion.id = 'abc'
		fusion.genes = (transcripts[junc1[0]].gene, transcripts[junc2[0]].gene)
		fusion.transcripts = (junc1[0], junc2[0])
		fusion.exons = (junc1[1][0][0] + 1, junc2[1][0][0] + 1)
		
		print 'kkko', fusion.as_tab()
		
		return fusion
	    
	return None
    
    def create_fusion_event(self, aligns, transcripts, junc1, junc2, pos=[], contig_breaks=None):
	if len(aligns) == 2:
	    fusion = call_event(aligns[0], aligns[1])
	else:
	    fusion = Adjacency((aligns[0].target, aligns[0].target),
	                       pos,
	                       '-',
	                       orients=('L', 'R'),
	                       contig_breaks = contig_breaks
	                       )
	                       
	fusion.rna_event = 'fusion'
	fusion.id = 'abc'
	fusion.genes = (transcripts[junc1[0]].gene, transcripts[junc2[0]].gene)
	fusion.transcripts = (junc1[0], junc2[0])
	fusion.exons = (junc1[1][0][0] + 1, junc2[1][0][0] + 1)
		
	return fusion
    
    def extract_novel_seq(self, adj):
	print 'ens', adj.contigs, adj.contig_breaks, self.contigs_fasta.fetch(adj.contigs[0])
	start, end = adj.contig_breaks if adj.contig_breaks[0] < adj.contig_breaks[1] else adj.contig_breaks[1], adj.contig_breaks[0]
	novel_seq = self.contigs_fasta.fetch(adj.contigs[0], start, end - 1)
	print 'ens', start, end, novel_seq, len(novel_seq)
	return novel_seq
	    
    def find_events(self, matches_by_transcript, align, transcripts):
	genes = Set([transcripts[txt].gene for txt in matches_by_transcript.keys()])
	#print 'genes', genes
	#print align.blocks
	#print align.query_blocks
	num_blocks = len(align.blocks)
	
	events = []
		
	print 'ttt', len(genes), genes
	if len(genes) > 1:
	    matches_by_block = {}
	    genes_in_block = {}
	    for i in range(num_blocks):
		matches_by_block[i] = {}
		genes_in_block[i] = Set()
		for transcript in matches_by_transcript.keys():
		    matches_by_block[i][transcript] = matches_by_transcript[transcript][i]
		    if matches_by_transcript[transcript][i] is not None:
			genes_in_block[i].add(transcripts[transcript].gene)
			
			
	    junctions = zip(range(num_blocks), range(num_blocks)[1:])
	    print 'tttk', junctions
	    for k in range(len(junctions)):
		i, j = junctions[k]
		print 'fff', i, j, junctions[k], genes_in_block[i], genes_in_block[j]
		if len(genes_in_block[i]) == 1 and len(genes_in_block[j]) == 1:
		    if list(genes_in_block[i])[0] != list(genes_in_block[j])[0]:
			print 'fusion splice', i, j		    
			junc1, junc2 = self.identify_fusion(matches_by_block[i], matches_by_block[j], transcripts)
			print 'fusion-genes', align.query, transcripts[junc1[0]].gene, transcripts[junc2[0]].gene
			pos = (align.blocks[i][1], align.blocks[j][0])
			contig_breaks = (align.query_blocks[i][1], align.query_blocks[j][0])
			event = self.create_fusion_event((align,), transcripts, junc1, junc2, pos=pos, contig_breaks=contig_breaks)
			events.append(event)
	    
		elif len(genes_in_block[i]) > 1:
		    print 'fusion block', i, j, genes_in_block
		    
		elif len(genes_in_block[j]) > 1:
		    print 'fusion block', i, j
		    
	    for i in genes_in_block.keys():
		if len(genes_in_block[i]) == 2:
		    print 'fusion block', i, matches_by_block[i]
		    self.identify_fusion_unknown_break(matches_by_block[i], transcripts)
	    		
	local_events = self.find_novel_junctions(matches_by_transcript, align, transcripts)
	if local_events:
	    for event in local_events:
		#event_type, gene, transcript, contigs=[], contig_blocks=[], contig_breaks=[], exons=[]
		print 'mmm', event
		adj = Adjacency((align.target,), event['pos'], '-', contig=align.query)
		if event['event'] in ('ins', 'del', 'dup', 'inv'):
		    adj.rearrangement = event['event']
		adj.rna_event = event['event']
		adj.genes = (list(genes)[0],)
		adj.transcripts = (event['transcript'][0],)
		adj.exons = event['exons'][0]
		adj.contig_breaks = event['contig_breaks']
		
		if adj.rearrangement == 'ins' or adj.rearrangement == 'dup':
		    novel_seq = self.extract_novel_seq(adj)
		    adj.novel_seq = novel_seq if align.strand == '+' else reverse_complement(novel_seq)
		    print 'ens2', adj.novel_seq
		    
		    if Event.is_ITD(adj, align, self.contigs_fasta):
			adj.rna_event = 'ITD'
			
		    adj.size = len(adj.novel_seq)
		    
		elif adj.rearrangement == 'del':
		    adj.size = adj.breaks[1] - adj.breaks[0] - 1
		
		events.append(adj)
		
		#events.append(Event(event['event'],
		                    #gene1 = list(genes)[0],
		                    #transcripts1 = event['transcript'],
		                    #gene2 = list(genes)[0],
		                    #transcripts2 = event['transcript'],
		                    #exons1 = event['exons'],
		                    #contigs = align.query,
		                    #contig_blocks = [event['blocks']],
		                    #pos1 = event['pos'][0],
		                    #pos2 = event['pos'][1],
		                    #chrom1 = align.target,
		                    #chrom2 = align.target
		                    #)
		              #)
	    
	return events
						    		    
	    
    def find_novel_junctions(self, block_matches, align, transcripts):	
	# find annotated junctions
	annotated = Set()
	for transcript in block_matches.keys():
	    matches = block_matches[transcript]
	    # sort multiple exon matches for single exon by exon num
	    [m.sort(key=lambda mm: int(mm[0])) for m in matches if m is not None]
	    
	    for i in range(len(matches) - 1):
		j = i + 1
		
		if matches[i] == None or matches[j] == None:
		    continue
		
		if (i, j) in annotated:
		    continue
		
		if self.is_junction_annotated(matches[i][-1], matches[j][0]):
		    annotated.add((i, j))
		    continue
		
	all_events = []
	for transcript in block_matches.keys():	  
	    matches = block_matches[transcript]
	    for i in range(len(matches) - 1):
		j = i + 1
			
		if matches[i] is None and matches[j] is not None:
		    # special case where the first 2 blocks is the utr and there's an insertion separating the 2 blocks
		    if i == 0:
			events = self.classify_novel_junction(matches[i], matches[j][0], chrom=align.target, blocks=align.blocks[i:j+1], transcript=transcripts[transcript])
			for e in events:
			    e['blocks'] = (i, j)
			    e['transcript'] = transcript
			all_events.extend(events)
		    continue
		
		# for retained intron, not a 'junction'
		if matches[i] is not None and len(matches[i]) > 1:
		    events = self.classify_novel_junction(matches[i], chrom=align.target, blocks=align.blocks[i], transcript=transcripts[transcript])
		    if events:
			for e in events:
			    e['blocks'] = (i, j)
			    e['transcript'] = transcript
			all_events.extend(events)
		    
		# skip if junction is annotated
		if (i, j) in annotated:
		    continue
		
		# for novel exon
		if matches[j] == None and matches[i] is not None:
		    for k in range(j + 1, len(matches)):
			if matches[k] is not None and matches[k][0][0] - matches[i][-1][0] == 1:
			    j = k
			    break
		
		if matches[i] is not None and matches[j] is not None:
		    # matches (i and j) are the flanking matches, blocks are the middle "novel" blocks
		    events = self.classify_novel_junction(matches[i][-1], matches[j][0], chrom=align.target, blocks=align.blocks[i:j+1], transcript=transcripts[transcript])
		    print 'event', i, j, align.blocks[i:j+1], transcript, events
		    if events:
			for e in events:
			    #e['blocks'] = (i, j)
			    e['blocks'] = range(i + 1, j)
			    print 'eee', e['blocks'], align.query_blocks
			    #e['contig_breaks'] = (align.query_blocks[i + 1][0], align.query_blocks[j][1])
			    e['contig_breaks'] = (align.query_blocks[i][1], align.query_blocks[j][0])
			    e['transcript'] = transcript
			    print 'eeef', e
			all_events.extend(events)
		
	if all_events:
	    grouped = {}
	    for event in all_events:
		key = str(event['blocks']) + event['event']
		try:
		    grouped[key].append(event)
		except:
		    grouped[key] = [event]
		
	    uniq_events = []
	    for events in grouped.values():
		#if len(events) == 1:
		    #uniq_events.append(events[0])
		#else:
		transcripts = [e['transcript'] for e in events]
		exons = [e['exons'] for e in events]
		
		contig_breaks = []
		if events[0].has_key('contig_breaks'):
		    contig_breaks = events[0]['contig_breaks']
		
		uniq_events.append({'event': events[0]['event'],
	                            'blocks': events[0]['blocks'],
	                            'transcript': transcripts,
	                            'exons': exons,
	                            'pos': events[0]['pos'],
	                            'contig_breaks': contig_breaks,
	                            }
	                           )
		    
	    print 'final', uniq_events
	    return uniq_events
	    
	
	
	
    def classify_novel_junction(self, match1, match2=None, chrom=None, blocks=None, transcript=None):
	events = []
	if type(blocks[0]) is int:
	    pos = (blocks[0], blocks[1])
	else:
	    pos = (blocks[0][1], blocks[1][0])
	print 'blocks', match1, match2, chrom, transcript, blocks, pos
	
	if match2 is None:
	    if len(match1) == 2:
		exons = [m[0] for m in match1]
		if match1[0][1] == '=>' and\
		   match1[-1][1] == '<=' and\
		   len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(match1) - 1:
		    events.append({'event': 'retained_intron', 'exons': exons, 'pos':pos})
		    
	elif match1 is None and match2 is not None:
	    if match2[0] == 0 and match2[1] == '<=':
		events.append({'event': 'ins', 'exons': match2[0], 'pos':pos})
		
	else:
	    if match2[0] > match1[0] + 1 and\
	       '=' in match1[1] and\
	       '=' in match2[1]:    
		events.append({'event': 'skipped_exon', 'exons': range(match1[0] + 1, match2[0]), 'pos':pos})
	       
	    if match1[0] == match2[0] and\
	       match1[1][1] == '<' and match2[1][0] == '>':
	       #match1[1] == '=<' and match2[1] == '>=':
		# examine gap: splice motif or not
		#print chrom, blocks
		#print blocks[0][1] + 1, blocks[0][1] + 2
		#print blocks[1][0] - 2, blocks[1][0] - 1
		#print self.ref_fasta.fetch(chrom, blocks[0][1], blocks[0][1] + 20)
		#print self.ref_fasta.fetch(chrom, blocks[1][0] - 3, blocks[1][0] + 27)
		#print transcript.strand
		
		expected_motif = 'gtag' if transcript.strand == '+' else 'ctac'
		motif = self.ref_fasta.fetch(chrom, blocks[0][1], blocks[0][1] + 2) + self.ref_fasta.fetch(chrom, blocks[1][0] - 3, blocks[1][0] - 1)
		#print blocks
		#print expected_motif, motif, blocks[1][0] - blocks[0][1] - 1
		
		gap_size = blocks[1][0] - blocks[0][1] - 1
		pos = (blocks[0][1], blocks[1][0])
		min_intron_size = 20
		event = None
		if gap_size > 0:
		    if gap_size > min_intron_size and expected_motif == motif.lower():
			event = 'novel_intron'
		    else:
			event = 'del'
		elif gap_size == 0:
		    event = 'ins'
		
		if event is not None:
		    events.append({'event': event, 'exons': [match1[0]], 'pos':pos})
		
	    if match1[1][1] == '=' and\
	       (match2[1][0] == '>=' or match2[1][0] == '<='):
		events.append({'event': 'alt_3', 'exons': [match1[0], match2[0]], 'pos':pos})
		
	    if (match1[1][1] == '>=' or match1[1][1] == '<=') and\
	       match2[1][0] == '=':
		events.append({'event': 'alt_5', 'exons': [match1[0], match2[0]], 'pos':pos})
		
	    if match2[0] == match1[0] + 1 and\
	       match1[1][1] == '=' and\
	       match2[1][0] == '=': 
		pos = ('%s-%s' % (blocks[0][1], blocks[1][0]), '%s-%s' % (blocks[-2][1], blocks[-1][0]))
		events.append({'event': 'novel_exon', 'exons': [], 'pos':pos})
	
	return events
		    
	
    def is_junction_annotated(self, match1, match2):
	if match2[0] == match1[0] + 1 and\
	   match1[1][1] == '=' and\
	   match2[1][0] == '=':
	    return True
	
	return False
		
		            
    def extract_aligns(self, alns):
	chimeric_aligns = {
	    'gmap': gmap.find_chimera,
	    }[self.aligner](alns, self.bam)
	if chimeric_aligns:
	    return chimeric_aligns
	else:
	    return [{
	        'gmap': gmap.find_single_unique,
	        }[self.aligner](alns, self.bam)]
        
    def match_exon(self, block, exon):
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
	
	#print block_matches[0][1][1] == '=', block_matches[0][1]
	#print block_matches[-1][1][0] == '='
	#print len([m for m in block_matches[1:-1] if m[1] == '==']) == len(block_matches) - 2
	#print len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(block_matches) - 1
	#print len([(a, b) for a, b in zip(exons, exons[1:]) if b == a - 1]) == len(block_matches) - 1	
		
	if block_matches[0][0][1][1] == '=' and\
	   block_matches[-1][0][1][0] == '=' and\
	   len([m for m in block_matches[1:-1] if m[0][1] == '==']) == len(block_matches) - 2 and\
	   (len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(block_matches) - 1 or\
	    len([(a, b) for a, b in zip(exons, exons[1:]) if b == a - 1]) == len(block_matches) - 1):
	    return True
		
	return False



def main(args, options):
    outdir = args[-1]
    # find events
    em = ExonMapper(*args)
    #events = em.match_all()
    Mapping.output(em.mappings, '%s/contig_mappings.tsv' % outdir)
    gene_mappings = Mapping.group(em.mappings)
    Mapping.output(gene_mappings, '%s/gene_mappings.tsv' % outdir)
    
    Event.output(em.events, outdir)
    
if __name__ == '__main__':
    usage = "Usage: %prog c2g_bam aligner contigs_fasta annotation_file genome_file(indexed) out_dir"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-b", "--r2c_bam", dest="r2c_bam_file", help="reads-to-contigs bam file")
    parser.add_option("-n", "--num_procs", dest="num_procs", help="Number of processes. Default: 5", default=5, type=int)
    parser.add_option("--junctions", dest="junctions", help="output junctions", action="store_true", default=False)
    
    (options, args) = parser.parse_args()
    if len(args) == 6:
        main(args, options)     