from shared.alignment import compare_chr, reverse_complement
from vcf import VCF
from sets import Set
import sys
import copy

class Variant:
    def __init__(self, event, adjs, chrom=None, pos=None):
	self.event = event
	self.adjs = adjs

	try:
	    self.id = '-'.join([str(adj.id) for adj in adjs])
	except:
	    self.id = None
	    
	self.chrom = chrom
	self.pos = pos
	self.filtered_out = False
	self.somatic = True
		
    def as_vcf(self, ref_fasta, genomic=True, size_threshold=100, insertion_as_sv=True):
	if len(self.adjs) == 2:
	    if self.event == 'INV':
		return self.inversion_as_vcf(ref_fasta)
	
	    elif self.event == 'TRL':
		return self.reciprocal_trl_as_vcf(ref_fasta)
	    
	    elif self.event == 'INS' and insertion_as_sv:
		return self.insertion_as_vcf(ref_fasta)
    
	output = []
	for adj in self.adjs:
	    output.append(adj.as_vcf(ref_fasta, size_threshold, genomic))
	return '\n'.join(output)
	        
    def inversion_as_vcf(self, ref_fasta):
	info = {}
	
	contigs = Set()
	for i in range(2):
	    for contig in self.adjs[i].contigs:
		contigs.add(contig)
	info['BKPTID'] = ','.join(contigs)
	
	if self.adjs[0].homol_seq and self.adjs[0].homol_seq[0]  != '-':
	    info['CIPOS'] = '0,%d' % len(self.adjs[0].homol_seq[0])
	
	if self.adjs[1].homol_seq and self.adjs[1].homol_seq[0]  != '-':
	    info['CIEND'] = '0,%d' % len(self.adjs[1].homol_seq[0])
	    
	return self.adjs[0].as_sv(ref_fasta, id_ext=self.id, info_ext=info)
        
    def reciprocal_trl_as_vcf(self, ref_fasta):
	parids = ((self.adjs[1].id + 'a', self.adjs[1].id + 'b'), 
	          (self.adjs[0].id + 'a', self.adjs[0].id + 'b'))
	
	event_id = self.event.upper() + self.id
	vcf_0 = self.adjs[0].as_breakends(ref_fasta, info_ext={'EVENT': (event_id, event_id),
	                                                       'PARID': parids[0]})
	vcf_1 = self.adjs[1].as_breakends(ref_fasta, info_ext={'EVENT': (event_id, event_id),
	                                                       'PARID': parids[1]})
	
	return vcf_0 + '\n' + vcf_1
    
    def insertion_as_vcf(self, ref_fasta):
	assert self.chrom is not None, 'must specify chrom for big insertion event'
	assert type(self.pos) is list or type(self.pos) is tuple, '"pos" for big insertion event must be list or tuple'
	info = {}
	
	contigs = Set()
	for i in range(2):
	    for contig in self.adjs[i].contigs:
		contigs.add(contig)
	info['BKPTID'] = ','.join(contigs)
	if self.adjs[0].insertion_size and self.adjs[0].insertion_size > 0:
	    info['SVLEN'] = self.adjs[0].insertion_size
	if self.pos[1] != self.pos[0]:
	    info['CIPOS'] = '0,%d' % (self.pos[1] - self.pos[0])
	    	
	return self.adjs[0].as_sv(ref_fasta, chrom_ext=self.chrom, pos_ext=self.pos[0], id_ext=self.id, info_ext=info)
    
    def check_adjs(self, adjs_passed):
	num_adjs_passed = 0
	for adj in self.adjs:
	    if adj.id in adjs_passed:
		num_adjs_passed += 1
		
	if num_adjs_passed == 2:
	    return True
	else:
	    return False


class Adjacency:
    def __init__(self, chroms, breaks, rearrangement, novel_seq='-',
                 contig=None, contig_breaks=None, contig_sizes=None, contig_support_span=None,
                 probes=None, orients=None, 
                 homol_seq=None, homol_coords=None,
                 aligns=None, align_types=None):
        """
        Args:
            chroms: (str, str) a tuple containing 2 strings (if list provided, will convert to tuple)
            breaks: (int, int) a tuple containing 2 ints (if list provided, will convert to tuple)
            rearrangement: (str) one of (del,ins,trl,inv,dup,trp)
            event: (str) one of (del,ins,trl,inv,trp,dup,ITD,PTD)
        """
        self.chroms = tuple(chroms)
        self.breaks = tuple(breaks)
        self.rearrangement = rearrangement
                        
        # for insertion event, or untemplated sequence at breakpoint
        self.novel_seq = novel_seq
        
        # for homologous sequences at breakpoint
        self.homol_seq = []
	if homol_seq is not None:
	    self.homol_seq.append(homol_seq)
	#else:
	    #self.homol_seq.append('-')
	self.homol_coords = []
	if homol_coords is not None:
	    self.homol_coords.append(homol_coords)
	#else:
	    #self.homol_coords.append(('-', '-'))
        
        # assume initialization always with single contig
        self.contigs = []
        if contig is not None:
            self.contigs.append(contig)
        self.contig_breaks = []
        if contig_breaks is not None:
            self.contig_breaks.append(contig_breaks)
	self.contig_support_span = []
        if contig_support_span is not None:
            self.contig_support_span.append(contig_support_span)
	self.contig_sizes = []
	if contig_sizes is not None:
            self.contig_sizes.append(contig_sizes)
	self.probes = []
	if probes is not None:
            self.probes.append(probes)
	self.aligns = []
	self.aligns.append(aligns)
	
	self.align_types = []
	self.align_types.append(align_types)
	
        # small-scale events may not have this information
        self.orients = orients
        	
	self.genes = ['NA', 'NA']
	self.transcripts = ['NA', 'NA']
	self.gene_strands = ['NA', 'NA']
	self.exons = ['NA', 'NA']
	self.exon_bounds = ['NA', 'NA']
	self.introns = ['NA', 'NA']
	self.rna_event = 'NA'
	self.fusion_type = 'NA'
	self.gene5 = 'NA'
	self.gene3 = 'NA'
	self.frames = 'NA'
			
	self.support = {'spanning':[], 'flanking':[], 'tiling':[]}
	self.support_total = {'flanking':'NA', 'spanning':'NA'}
	self.support_final = None
	
	self.support_normal = {'spanning':[], 'flanking':[], 'tiling':[]}
	self.support_total_normal = {'flanking':'NA', 'spanning':'NA'}
	self.support_final_normal = None
	
	self.filtered_out = None
	self.filtered_out_normal = None
	
	self.somatic = None
	self.insertion_size = None
	
	self.partner_contig = None
	
	self.stigmas = Set()
	
	self.dubious = False
		
    def debug(self):
        print '%s %s %s %s %s %s' % (self.rearrangement, self.chroms, self.breaks, self.get_size(), ','.join(self.contigs), self.orients)
        
    def get_size(self):
        size = 0
        if self.align_types[0] == 'gapped' and self.rearrangement in ('ins', 'dup'):
            if self.novel_seq:
                size = len(self.novel_seq)
	elif self.insertion_size is not None:
	    size = self.insertion_size
	elif self.rearrangement == 'ins':
	    if self.novel_seq != 'NA':
		size = len(self.novel_seq)
	elif self.rearrangement == 'del':
	    size = self.breaks[1] - self.breaks[0] - 1
        elif self.rearrangement != 'trl':
            size = self.breaks[1] - self.breaks[0] + 1
	else:
	    size = 'NA'
        
        return size
        
    def as_vcf(self, ref_fasta, size_threshold, genomic=True):	    
	size = self.get_size()
	if (self.rearrangement == 'ins' or self.rearrangement == 'del') and type(size) is int and size <= size_threshold:
	    return self.as_indel(ref_fasta)
	elif self.rearrangement == 'trl':
	    return self.as_breakends(ref_fasta, genomic)
	else:
	    return self.as_sv(ref_fasta)
        
    def as_breakends(self, ref_fasta, genomic=True, max_novel_seq_len=50, info_ext=None, parids=None, event=None):
        chroms = map(lambda c: c.lstrip('chr'), self.chroms)
	alt_chroms = chroms[:]
	pos = list(self.breaks)
	alt_pos = pos[:]
	# inserted novel sequences
        inserted_seqs = ['','']
	if self.novel_seq and self.novel_seq != 'NA' and self.novel_seq != '-':
	    if len(self.novel_seq) > max_novel_seq_len:
		alt_chroms[0] = '<%s>' % self.contigs[0]
		alt_chroms[1] = '<%s>' % self.contigs[0]
		alt_pos[1] = self.contig_breaks[0][0] + 1
		alt_pos[0] = self.contig_breaks[0][1] - 1
	    else:
		if len(self.aligns[0]) == 1:
		    inserted_seqs[0] = self.novel_seq if self.aligns[0][0].strand == '+' else reverse_complement(self.novel_seq)
		    inserted_seqs[1] = self.novel_seq if self.aligns[0][0].strand == '+' else reverse_complement(self.novel_seq)
		else:
		    inserted_seqs[0] = self.novel_seq if self.aligns[0][0].strand == '+' else reverse_complement(self.novel_seq)
		    inserted_seqs[1] = self.novel_seq if self.aligns[0][1].strand == '+' else reverse_complement(self.novel_seq)
		
	# microhomology, cipos
	cipos = None
	homol_len = None
	homol_seq = None
	if self.homol_seq and self.homol_seq[0] != '-' and len(self.homol_seq) > 0:
	    homol_seq = self.homol_seq[0].upper()
	    homol_len = len(self.homol_seq[0])
	    contig_breaks = self.contig_breaks[0]
	    # e.g. GMAP
	    if contig_breaks[0] + 1 == contig_breaks[1]:
		pass
	    # e.g. BWA-mem
	    elif contig_breaks[0] >= contig_breaks[1]:
		pos[0] -= homol_len
		alt_pos[1] += homol_len
		cipos = '0,%d' % homol_len
	    
        refs = (ref_fasta.fetch(self.chroms[0], self.breaks[0] - 1, self.breaks[0]).upper(),
                ref_fasta.fetch(self.chroms[1], self.breaks[1] - 1, self.breaks[1]).upper())

        ids = ('%s%s' % (self.id, 'a'),
               '%s%s' % (self.id, 'b'))
	        
        svtype = 'BND' if genomic else 'FND'
        infos = [{'SVTYPE':svtype, 'MATEID':ids[1], 'EVENTTYPE':self.rearrangement.upper()},
                 {'SVTYPE':svtype, 'MATEID':ids[0], 'EVENTTYPE':self.rearrangement.upper()}]
	if cipos is not None:
	    infos[0]['CIPOS'] = cipos
	    infos[1]['CIPOS'] = cipos
	if homol_len is not None:
	    infos[0]['HOMLEN'] = homol_len
	    infos[1]['HOMLEN'] = homol_len
	if homol_seq is not None:
	    infos[0]['HOMSEQ'] = homol_seq
	    infos[1]['HOMSEQ'] = homol_seq
	    
	# read support
	if self.support_final is not None:
	    infos[0]['READSUPPORT'] = self.support_final
	    infos[1]['READSUPPORT'] = self.support_final
	    
	adj_size = self.get_size()
	if type(adj_size) is int:
	    infos[0]['SVLEN'] = adj_size
	    infos[1]['SVLEN'] = adj_size
	    
	# somatic
	if self.somatic:
	    infos[0]['SOMATIC'] = 'SOMATIC'
	    infos[1]['SOMATIC'] = 'SOMATIC'
	    
        # contig and contig breakpoints
        if self.contigs:
            for i in range(2):
                infos[i]['BKPTID'] = ','.join(self.contigs)
            
        if self.contig_breaks and len(self.contig_breaks) == len(self.contigs):
            contig_breaks = []
            for bk in self.contig_breaks:
                if len(bk) == 2:
                    contig_breaks.append('%s-%s' % (bk[0], bk[1]))
                else:
                    print 'error'
                    
            if len(contig_breaks) == len(self.contigs):
                for i in range(2):
                    infos[i]['CTG_BKS'] = ','.join(contig_breaks)
		    	    
	# external info - overrides given info
	if info_ext:
	    for key, value in info_ext.iteritems():
		if len(value) == 2:
		    infos[0][key] = value[0]
		    infos[1][key] = value[1]
        
        if self.orients[0] == 'L':
            # LL
            if self.orients[1] == 'L':
                alts = ('%s%s]%s:%s]' % (refs[0], inserted_seqs[0], alt_chroms[1], alt_pos[1]),
                        '%s%s]%s:%s]' % (refs[1], inserted_seqs[1], alt_chroms[0], alt_pos[0]))
            # LR
            else:
                alts = ('%s%s[%s:%s[' % (refs[0], inserted_seqs[0], alt_chroms[1], alt_pos[1]),
                        ']%s:%s]%s%s' % (alt_chroms[0], alt_pos[0], inserted_seqs[1], refs[1]))
        else:
            # RL
            if self.orients[1] == 'L':
                alts = (']%s:%s]%s%s' % (chroms[1], alt_pos[1], inserted_seqs[0], refs[0]),
                        '%s%s[%s:%s[' % (refs[1], inserted_seqs[1], chroms[0], alt_pos[0])) 
            # RR
            else:
                alts = ('[%s:%s[%s%s' % (chroms[1], alt_pos[1], inserted_seqs[0], refs[0]),
                        '[%s:%s[%s%s' % (chroms[0], alt_pos[0], inserted_seqs[1], refs[1]))
	
        breakends = map(lambda i: '\t'.join([chroms[i], str(pos[i]), ids[i], refs[i], alts[i], '.', '.', VCF.info_dict_to_str(infos[i])]), range(2))
        
        return '\n'.join(breakends)
                
    def as_indel(self, ref_fasta):
        chrom = self.chroms[0].lstrip('chr')
	pos = self.breaks[0]
	
	ref = alt = None
	size = self.get_size()
	if self.rearrangement == 'del':
	    ref = ref_fasta.fetch(self.chroms[0], self.breaks[0] - 1, self.breaks[1] - 1).upper()
	    alt = ref_fasta.fetch(self.chroms[0], self.breaks[0] - 1, self.breaks[0]).upper()
	    
	else:
	    ref = ref_fasta.fetch(self.chroms[0], self.breaks[0] - 1, self.breaks[1]).upper()
	    alt = ref + self.novel_seq.upper()
	    
	id = self.id
	qual = '.'
	filter = '.'
	info = {
	        'BKPTID':','.join(self.contigs),
	        }
	
	# read support
	if self.support_final is not None:
	    info['READSUPPORT'] = self.support_final
	    
	# somatic
	if self.somatic:
	    info['SOMATIC'] = 'SOMATIC'
	    
	if ref is not None and alt is not None:
	    fields = [chrom, pos, id, ref, alt, qual, filter, VCF.info_dict_to_str(info)]
	    return '\t'.join(map(str, fields))
	
    def as_sv(self, ref_fasta, id_ext=None, info_ext=None, chrom_ext=None, pos_ext=None):
        chrom = self.chroms[0] if chrom_ext is None else chrom_ext
	pos = self.breaks[0] if pos_ext is None else pos_ext
	
	chrom = chrom.lstrip('chr')
	
	alt = None
	ref = ref_fasta.fetch(self.chroms[0], self.breaks[0] - 1, self.breaks[0]).upper()
	sv_len = self.get_size()
	end = None
	if type(sv_len) is int:
	    end = pos + sv_len
	    
	if self.rearrangement == 'del':
	    alt = '<DEL>'
	    sv_type = 'DEL'
	    if type(sv_len) is int:
		sv_len = -1 * sv_len
		end = pos - sv_len
	    
	elif self.rearrangement == 'dup':
	    alt = '<DUP:TANDEM>'
	    sv_type = 'DUP'
	    
	elif self.rearrangement == 'inv':
	    alt = '<INV>'
	    sv_type = 'INV'
	    
	elif self.rearrangement == 'ins':
	    alt = '<INS>'
	    sv_type = 'INS'
	    end = pos
	    
	id = self.id if id_ext is None else id_ext
	qual = '.'
	filter = '.'
	info = {'SVTYPE': sv_type,
	        'END': end,
	        'BKPTID':','.join(self.contigs),
	        }
	if end is not None:
	    info['END'] = end
	if type(sv_len) is int:
	    info['SVLEN'] = sv_len
	
	# read support
	if self.support_final is not None:
	    info['READSUPPORT'] = self.support_final
	    
	# somatic
	if self.somatic:
	    info['SOMATIC'] = 'SOMATIC'
    
	cipos = None
	homol_len = None
	homol_seq = None
	if self.homol_seq and self.homol_seq[0] != '-':
	    homol_seq = self.homol_seq[0].upper()
	    homol_len = len(self.homol_seq[0])
	    contig_breaks = self.contig_breaks[0]
	    # e.g. GMAP
	    if contig_breaks[0] + 1 == contig_breaks[1]:
		#print 'gmap', contig_breaks
		pass
	    # e.g. BWA-mem
	    elif contig_breaks[0] >= contig_breaks[1]:
		cipos = '0,%d' % homol_len
		
	if cipos is not None:
	    info['CIPOS'] = cipos
	    info['CIPOS'] = cipos
	if homol_len is not None:
	    info['HOMLEN'] = homol_len
	    info['HOMLEN'] = homol_len
	if homol_seq is not None:
	    info['HOMSEQ'] = homol_seq
	    info['HOMSEQ'] = homol_seq
	    
	# external info - overrides given info
	if info_ext:
	    for key, value in info_ext.iteritems():
		if key == 'SVLEN' and value == 'NA':
		    continue
		info[key] = value
	
	if ref is not None and alt is not None:
	    fields = [chrom, pos, id, ref, alt, qual, filter, VCF.info_dict_to_str(info)]
	    return '\t'.join(map(str, fields))
        
    @classmethod
    def show_tab_headers(cls):
        headers = ('ID',
                   'event',
                   'chrom1',
                   'break1',
                   'orient1',
                   'chrom2',
                   'break2',
                   'orient2',
                   'size',
                   'contig',
                   'contig-break1',
                   'contig-break2',
	           'homol_seq',
	           'homol_coord1',
	           'homol_coord2',
	           'novel_seq',
	           'probe',
	           'gene1',
	           'gene2',
	           'transcript1',
	           'transcript2',
	           'exon1',
	           'exon2',
	           'exon_bound1',
	           'exon_bound2',
	           'intron1',
	           'intron2',
	           'RNA_event',
	           'fusion_type',
	           '5_gene',
	           '3_gene',
	           'frame',
	           'spanning_reads',
	           'flanking_pairs',
	           'read_tiling',
                   )
        return '\t'.join(headers)
    
    def as_tab(self, expand_contigs=False):        
        outputs = []
	
	if not expand_contigs:
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
	    data.append(','.join(self.homol_seq))
	    try:
		data.append(','.join([str(b[0]) for b in self.homol_coords]))
	    except:
		data.append('-')
	    try:
		data.append(','.join([str(b[1]) for b in self.homol_coords]))
	    except:
		data.append('-')
	    data.append(self.novel_seq)
	    try:
		data.append(self.probes[0])
	    except:
		data.append('-')
	    data.append(self.genes[0])
	    data.append(self.genes[1])
	    data.append(self.transcripts[0])
	    data.append(self.transcripts[1])
	    data.append(str(self.exons[0]))
	    data.append(str(self.exons[1]))
	    data.append(str(self.exon_bounds[0]))
	    data.append(str(self.exon_bounds[1]))
	    data.append(str(self.introns[0]))
	    data.append(str(self.introns[1]))
	    data.append(self.rna_event)
	    data.append(self.fusion_type)
	    data.append(self.gene5)
	    data.append(self.gene3)
	    data.append(str(self.frames))
	    data.append(str(self.support_total['spanning']))
	    data.append(str(self.support_total['flanking']))
	    data.append(','.join(str(t) for t in self.support['tiling']))
	    outputs.append('\t'.join(map(str, data)))
	    
	else:
	    for i in range(len(self.contigs)):
		data = []
		data.append(self.id)
		data.append(self.rearrangement)
		for j in range(2):
		    data.append(self.chroms[j])
		    data.append(str(self.breaks[j]))
		    data.append(self.orients[j])
		data.append(str(self.get_size()))
		data.append(self.contigs[i])
		data.append(str(self.contig_breaks[i][0]))
		data.append(str(self.contig_breaks[i][1]))
		data.append(str(self.homol_seqs[i]))
		data.append(str(self.homol_coords[i][0]))
		data.append(str(self.homol_coords[i][1]))
		data.append(self.novel_seq)  
		data.append(self.probes[i])
		data.append(self.genes[0])
		data.append(self.genes[1])
		data.append(self.transcripts[0])
		data.append(self.transcripts[1])
		data.append(str(self.exons[0]))
		data.append(str(self.exons[1]))
		data.append(str(self.introns[0]))
		data.append(str(self.introns[1]))
		data.append(str(self.exon_bounds[0]))
		data.append(str(self.exon_bounds[1]))
		data.append(self.rna_event)
		data.append(self.fusion_type)
		data.append(self.gene5)
		data.append(self.gene3)
		data.append(str(self.frames))
		if len(self.read_supports) == len(self.contigs):
		    data.append(str(self.support['spanning'][i]))
		    data.append(str(self.support['flanking'][i]))
		    data.append(str(self.support['tiling'][i]))
		outputs.append('\t'.join(map(str, data)))
            
        return '\n'.join(outputs)
    
    def as_bed(self):
	outputs = []
	if len(self.chroms) == 2:
	    for i in range(len(self.chroms)):
		data = []
		data.append(self.chroms[i])
		data.append(str(self.breaks[i] - 1))
		data.append(str(self.breaks[i]))
		data.append(self.key())
		data.append('0')
		data.append('+')
		outputs.append('\t'.join(data))
		
	else:
	    data = []
	    data.append(chrom)
	    data.append(str(self.breaks[0] - 1))
	    data.append(str(self.breaks[1]))
	    data.append(self.key())
	    data.append('0')
	    data.append('+')
	    outputs.append('\t'.join(data))
            
	return '\n'.join(outputs)
    
    def as_bedpe(self):
	data = []
	for i in range(len(self.chroms)):
	    data.append(self.chroms[i])
	    data.append(str(self.breaks[i] - 1))
	    data.append(str(self.breaks[i]))
	data.append(self.key())
	data.append('.')
	for i in range(len(self.orients)):
	    if i == 'L':
		data.append('+')
	    else:
		data.append('-')
		
	return '\t'.join(data)
    
    def contig_support_span(self, contig_index):
	if self.homol_coords and self.homol_coords[contig_index][0] is int and self.homol_coords[contig_index][1] is int:
	    return (self.homol_coords[contig_index][0], self.homol_coords[contig_index][1])
	else:
	    return (self.contig_breaks[contig_index][0], self.contig_breaks[contig_index][1])
    
	    
    def sum_support(self, normal=False):
	(support, support_total) = (self.support, self.support_total) if not normal else (self.support_normal, self.support_total_normal)
	for kind, nums in support.iteritems():
	    if kind == 'tiling':
		continue
	    
	    support_total[kind] = sum(nums)
	            
    def key(self, transcriptome=False, include_novel_seq=False):
	"""Constructs a unique key for grouping adjacencies
	
	Args:
	    transcriptome: (boolean) whether adjacency is genomic or transcriptomic
	Returns:
	    A string that is used for grouping adjacencies
	"""
	info = [self.rearrangement]
	if transcriptome:
	    info = [self.rna_event]
	for i in (0,1):
	    info.append(self.chroms[i])
	    info.append(self.breaks[i])
	    if self.orients:
		info.append(self.orients[i])
	    if include_novel_seq:
		info.append(self.novel_seq)
	return '-'.join(map(str, info))
	        
    @classmethod
    def extract_probe(cls, contig_seq, contig_breaks, len_on_each_side=25):
	start, end = contig_breaks
	if contig_breaks[0] > contig_breaks[1]:
	    start, end = contig_breaks[1], contig_breaks[0]
		
	upstream_coord = max(0, start - len_on_each_side)
	downstream_coord = min(end - 1 + len_on_each_side, len(contig_seq))
	    
	return contig_seq[upstream_coord:downstream_coord], start - upstream_coord
    
    @classmethod
    def extract_probe_new(cls, contig_seq, contig_breaks, len_on_each_side=50, kmer_size=None, min_buffer=1):
	probe = 'NA'
	contig_breaks_sorted = sorted(contig_breaks)
	print contig_breaks_sorted, min_buffer, kmer_size	
	start = contig_breaks_sorted[1] + min_buffer - kmer_size + 1
	end = contig_breaks_sorted[0] - min_buffer + kmer_size - 1
	
	probe_size = 0
	if start >= 1 and end <= len(contig_seq) and end - start + 1 >= kmer_size:
	    probe = contig_seq[start - 1 : end]
	    probe_size = len(probe)
	    	
	return probe
	
    def extract_subseqs(self, contig_fasta):
	subseqs = []
	aligns = self.aligns[0]
	for i in range(len(aligns)):
	    subseqs.append(contig_fasta.fetch(self.contigs[0], aligns[i].qstart - 1, aligns[i].qend))
	
	return subseqs
    
    @classmethod
    def group_inversions(cls, adjs):
	"""Group 2 inversion adjacencies into a single event"""
	inversions = sorted(adjs, key=lambda adj: (adj.chroms[0], adj.breaks[0]))
	
	max_homology = 25
	variants = []
	i = 0
	while i < len(inversions) - 1:
	    if inversions[i].chroms[0] == inversions[i + 1].chroms[0] and\
	       inversions[i + 1].breaks[0] - inversions[i].breaks[0] <= max_homology and\
	       ((inversions[i].orients == ('L', 'L') and inversions[i + 1].orients == ('R', 'R')) or
	        (inversions[i].orients == ('R', 'R') and inversions[i + 1].orients == ('L', 'L'))):
		
		(adj1, adj2) = (inversions[i], inversions[i + 1]) if inversions[i].orients == ('L', 'L') else (inversions[i + 1], inversions[i])
		
		variants.append(Variant('INV', [adj1, adj2]))
		i += 2
		
	    else:
		if not inversions[i].dubious:
		    variants.append(Variant('INV', [inversions[i]]))
		i += 1
		
	if i == len(inversions) - 1 and not inversions[i].dubious:
	    variants.append(Variant('INV', [inversions[i]]))
	    
	return variants
    
    @classmethod
    def extract_interchrom_ins(cls, trls):
	"""Group translocation events that are possibly insertions
	
	Args:
	    trls: (List) Adjacencies that are translocations
	Returns:
	    A list of variants that are possibly insertions
	    A list of remaining translocation adjacencies
        """
	def pair_up(trls, ins_at_first=True):
	    """Pairs up 2 translocation adjacencies that can be potentially insertions
	    
	    Conditions for having an insertion event:
	    1. At most only 2 chromosomes are invovled: source and target
	    2. the break positions on the target chromosome is reasonably close (500bp)
	    3. the orientations and positions of the breaks on the source chromosome makes sense
	       i.e. if it is (L,R), then the first break pos must be greater than the second
		    if it is (R,L), then the second break pos must be greater than the first
	    
	    Args:
	        trls: (List) Adjacencies that are translocations
	        ins_at_first: (Boolean) Assuming first break position is target (second break is source)
	    Returns:
	        A list of variants that are possibly insertions
		A list of remaining translocation adjacencies
	    """
	    neighborhood= 500
	    insertions = []
	    used = Set()
	    for i in range(len(trls)):
		if i in used:
		    continue
		
		for j in range(len(trls)):
		    if i == j or j in used:
			continue
		    
		    if trls[i].chroms[0] != trls[j].chroms[0] or trls[i].chroms[1] != trls[j].chroms[1]:
			continue
		    
		    anchor_dubious = False
		    if ins_at_first:
			target_chrom = trls[i].chroms[0]
			target_breaks = trls[i].breaks[0], trls[j].breaks[0]
			target_orients = trls[i].orients[0], trls[j].orients[0]
			source_chroms = trls[i].chroms[1]
			source_breaks = trls[i].breaks[1], trls[j].breaks[1]
			source_orients = trls[i].orients[1], trls[j].orients[1]
			
			if trls[i].aligns[0][0].dubious:
			    print 'anchor', target_chrom, target_breaks, trls[i].aligns[0][0].target, trls[i].aligns[0][0].tstart, trls[i].aligns[0][0].tend
			    anchor_dubious = True
		    else:
			target_chrom = trls[i].chroms[1]
			target_breaks = trls[i].breaks[1], trls[j].breaks[1]
			target_orients = trls[i].orients[1], trls[j].orients[1]
			source_chroms = trls[i].chroms[0]
			source_breaks = trls[i].breaks[0], trls[j].breaks[0]
			source_orients = trls[i].orients[0], trls[j].orients[0]
			
			if trls[i].aligns[0][1].dubious:
			    #print 'anchor', target_chrom, target_breaks, trls[i].aligns[0][1].target, trls[i].aligns[0][1].tstart, trls[i].aligns[0][1].tend
			    anchor_dubious = True
			    
		    #print 'anchor', anchor_dubious
		    if anchor_dubious:
			continue
								    
		    if abs(target_breaks[0] - target_breaks[1]) <= neighborhood and\
		       target_orients[0] != target_orients[1] and\
		       ((source_orients[0] == 'L' and source_orients[1] == 'R' and source_breaks[0] > source_breaks[1]) or\
		       (source_orients[0] == 'R' and source_orients[1] == 'L' and source_breaks[1] < source_breaks[0])):
			
			if source_breaks[0] > source_breaks[1]:
			    insertion_size = source_breaks[0] - source_breaks[1] + 1
			else:
			    insertion_size = source_breaks[1] - source_breaks[0] + 1
			
			trls[i].rearrangement = 'ins'		    
			trls[j].rearrangement = 'ins'
			trls[i].insertion_size = insertion_size
			trls[j].insertion_size = insertion_size
			
			insertion = Variant('INS', [trls[i], trls[j]], chrom=target_chrom, pos=target_breaks)
			insertions.append(insertion)
			used.add(i)
			used.add(j)
						
	    unused = [trls[i] for i in range(len(trls)) if not i in used]
	    return insertions, unused
	
	variants = []
	
	# if insertion is at first chromosome
	insertions, trls_remained = pair_up(trls)
	variants.extend(insertions)
		
	# if insertion is at second chromosome
	if trls_remained:
	    insertions, trls_remained = pair_up(trls_remained, ins_at_first=False)
	    variants.extend(insertions)
			
	return variants, trls_remained
    
    @classmethod
    def group_trls(cls, adjs):
	"""Group 2 translocation adjacencies into single reciprocal event"""
	trls = sorted([adj for adj in adjs if not adj.dubious], key=lambda adj: (adj.chroms[0], adj.breaks[0]))
	
	neighborhood = 10000
	variants = []
	i = 0
	if len(trls) > 1:
	    while i <= len(trls) - 1:
		if i == len(trls) - 1:
		    variants.append(Variant('TRL', [trls[i]]))
		    break
		
		elif trls[i].chroms[0] == trls[i + 1].chroms[0] and\
		   trls[i].chroms[1] == trls[i + 1].chroms[1] and\
		   abs(trls[i + 1].breaks[0] - trls[i].breaks[0]) <= neighborhood and\
		   abs(trls[i + 1].breaks[1] - trls[i].breaks[1]) <= neighborhood and\
		   ((trls[i].orients == ('L', 'R') and trls[i + 1].orients == ('R', 'L')) or\
		    (trls[i].orients == ('R', 'L') and trls[i + 1].orients == ('L', 'R')) or\
		    (trls[i].orients == ('L', 'L') and trls[i + 1].orients == ('R', 'R')) or\
		    (trls[i].orients == ('R', 'R') and trls[i + 1].orients == ('L', 'L'))
		    ):
		    
		    variants.append(Variant('TRL', [trls[i], trls[i + 1]]))
		    i += 2
		    
		else:
		    if not trls[i].dubious:
			variants.append(Variant('TRL', [trls[i]]))
		    i += 1
		    
	elif len(trls) == 1 and not trls[0].dubious:
	    variants.append(Variant('TRL', [trls[0]]))
		
	return variants
    
    @classmethod
    def extract_imprecise_ins(cls, adjs, debug=False):
	def screen_insertion(adj1, adj2):
	    neighborhood= 50
	    
	    # s1, t1 = source can 0, 1 of adj1.breaks
	    for s1 in range(2):
		t1 = 1 if s1 == 0 else 0
		for s2 in range(2):
		    t2 = 1 if s2 == 0 else 0
		    
		    if adj1.chroms[t1] == adj2.chroms[t2] and\
		       abs(adj1.breaks[t1] - adj2.breaks[t2]) <= neighborhood and\
		       adj1.orients[t1] != adj2.orients[t2]:
			if adj1.chroms[s1] == adj2.chroms[s2]:
			    if abs(adj1.breaks[s1] - adj2.breaks[s2]) > neighborhood and\
			       adj1.orients[s1] != adj2.orients[s2] and\
			       adj1.orients[s1] != adj1.orients[t1] and\
			       adj2.orients[s2] != adj2.orients[t2]:
				return (s1, s2, t1, t2)
			    
			# source coming from different chromosomes or other genomes
			# can only just check the target	
			else:
			    return (s1, s2, t1, t2)		
			    
	    return None
	    
	insertions = []
	used_contigs = Set()
	used_adjs = Set()
	for i in range(len(adjs) - 1):	
	    for j in range(i + 1, len(adjs)):		
		# same contig - skip
		if adjs[i].contigs[0] == adjs[j].contigs[0]:
		    continue
		
		if adjs[i].contigs[0] in used_contigs or adjs[j].contigs[0] in used_contigs:
		    continue
		
		ins_event = screen_insertion(adjs[i], adjs[j])
		    
		if ins_event is not None:
		    src_chroms = adjs[i].chroms[ins_event[0]], adjs[j].chroms[ins_event[1]]
		    src_coords = sorted([adjs[i].breaks[ins_event[0]], adjs[j].breaks[ins_event[0]]])
		    target_chrom = adjs[i].chroms[ins_event[2]]
		    target_coords = sorted([adjs[i].breaks[ins_event[2]], adjs[j].breaks[ins_event[3]]])
		    ins_size = None
		    if src_chroms[0] == src_chroms[1]:
			ins_size = src_coords[1] - src_coords[0] + 1
		    if debug:
			sys.stdout.write('ins_grouped contigs:%s,%s source:%s:%s %s:%s target:%s:%s-%s\n' % (adjs[i].contigs[0], 
			                                                                                  adjs[j].contigs[0], 
			                                                                                  src_chroms[0],
			                                                                                  src_coords[0],
			                                                                                  src_chroms[1],
			                                                                                  src_coords[1],
			                                                                                  target_chrom,
			                                                                                  target_coords[0],
			                                                                                  target_coords[1])
			                 )
		    adjs[i].rearrangement = 'ins'
		    adjs[j].rearrangement = 'ins'
		    adjs[i].insertion_size = ins_size
		    adjs[j].insertion_size = ins_size
		    
		    insertions.append(Variant('INS', (adjs[i], adjs[j]), chrom=target_chrom, pos=target_coords))
		    
		    used_contigs.add(adjs[i].contigs[0])
		    used_contigs.add(adjs[j].contigs[0])
		    
		    used_adjs.add(adjs[i].id)
		    used_adjs.add(adjs[j].id)
		
	return insertions, used_adjs    	
		        
    @classmethod
    def merge(cls, adjs, transcriptome=False):
	"""Merge adjacencies that have the same breakpoint (and same event type) together
	Args:
	    adjs: (list) Adjacency
	    transcriptome: (boolean) whether adjacency is genomic or transcriptomic
	Returns:
	    List of adjs with subsets that represent the same adjacency merged
	"""
	keys = {}
	for adj in adjs:
	    key = adj.key(transcriptome=transcriptome)

	    if not keys.has_key(key):
		keys[key] = copy.deepcopy(adj)
	    else:
		first_adj = keys[key]
		first_adj.contigs.append(adj.contigs[0])
		first_adj.contig_breaks.append(adj.contig_breaks[0])
		first_adj.contig_sizes.append(adj.contig_sizes[0])
		first_adj.contig_support_span.append(adj.contig_support_span[0])
		if adj.probes:
		    first_adj.probes.append(adj.probes[0])
		else:
		    first_adj.probes.append('-')
		first_adj.aligns.append(adj.aligns[0])
		first_adj.align_types.append(adj.align_types[0])
		if adj.homol_seq:
		    first_adj.homol_seq.append(adj.homol_seq[0])
		if adj.homol_coords:
		    first_adj.homol_coords.append(adj.homol_coords[0])
		for support_type in ('spanning', 'tiling', 'flanking'):
		    if adj.support[support_type]:
			first_adj.support[support_type].append(adj.support[support_type][0])
	
	# for generating ID
	count = 1
	adjs_merged = []
	for key, adj in sorted(keys.iteritems()):
	    adj.id = str(count)
	    #if len(adj.contigs) > 1:
		#if adj.homol_seq and adj.homol_coords:
		    #adj.contig_sizes, adj.contigs, adj.contig_breaks, adj.probes, adj.aligns, adj.align_types, adj.homol_seq, adj.homol_coords, adj.support['spanning'] =\
		           #(list(t) for t in zip(*sorted(zip(adj.contig_sizes, 
		                                             #adj.contigs, 
		                                             #adj.contig_breaks, 
		                                             #adj.probes,
		                                             #adj.aligns,
		                                             #adj.align_types,
		                                             #adj.homol_seq,
		                                             #adj.homol_coords,
		                                             #adj.support['spanning']), 
		                                         #reverse=True)))
		#else:
		    #adj.contig_sizes, adj.contigs, adj.contig_breaks, adj.probes, adj.aligns, adj.align_types, adj.support['spanning'] =\
		           #(list(t) for t in zip(*sorted(zip(adj.contig_sizes, 
		                                             #adj.contigs, 
		                                             #adj.contig_breaks, 
		                                             #adj.probes,
		                                             #adj.aligns,
		                                             #adj.align_types,
		                                             #adj.support['spanning']),
		                                         #reverse=True)))	    	
	    count += 1
	    adjs_merged.append(adj)
	    		
	return adjs_merged
        
    @classmethod
    def realign(cls, adjs, out_dir, aligner, 
                probe=False, subseq=False, 
                contigs_fasta = None,
                use_realigns=False, 
                name_sep='-', genome=None, index_dir=None, num_procs=None):
	"""Aligns probe and subsequences against reference genome
	
	The output of the fasta sequences will be called 'realign.fa',
	and the alignments will be in 'realign.bam', put in the output directory
	
	Args:
	    adjs: (list) Adjacencies for extracting sequences
	    out_dir: (str) full path of output directory for storing sequences and alignments
	    aligner: (str) Name of aligner (gmap, bwa_mem)
	    probe: (boolean) Align probe sequence. Default: None
	    subseq: (boolean) Align subsequences. Default: None
	    contigs_fasta: (pysam.fastafile) For extracting sub-sequences
	    use_realigns: (boolean) Use existing realignments. Default: False
	    name_sep: (str) Character used to combine various info into query name
	    genome: (str) prefix of the index of the refernece genome
	    index_dir: (str) full path of the directory of location of the genome index
	    num_procs: (int) number of threads to run the alignment
	
	Returns:
	    pysam.samfile handle of realignment bam
	"""
	def write_probe(adj, out, name_sep):
	    """Outputs the probe sequence to output file
	    
	    Args:
		adj: Variant object
		out: Filehandle of output file
		name_sep: (str) Character used to combine various info into query name
	    """
	    out.write('>%s%s%s\n%s\n' % (adj.contigs[0], name_sep, adj.key(), adj.probes[0]))
	    
	def write_subseq(adj, out, name_sep, contigs_fasta):
	    """Outputs the sub-sequence of a split alignment to output file
	
	    Args:
	        adj: Variant object
		out: Filehandle of output file
		name_sep: (str) Character used to combine various info into query name
	    """	    
	    subseqs = adj.extract_subseqs(contigs_fasta)
	    for i in range(len(subseqs)):
		out.write('>%s%s%s%s%d\n%s\n' % (adj.contigs[0], name_sep, adj.key(), name_sep, i, subseqs[i]))
	    
	from shared import gmap
	from shared import bwa_mem
	
	prefix = 'realign'
	if not use_realigns:
	    out_file = '%s/%s.fa' % (out_dir, prefix)
	    out = open(out_file, 'w')
	    for adj in adjs:
		if probe:
		    write_probe(adj, out, name_sep)
		if subseq and contigs_fasta is not None:
		    write_subseq(adj, out, name_sep, contigs_fasta)
	    out.close()
	    
	# run aligner
	realign_bam_file = '%s/%s.bam' % (out_dir, prefix)
	if not use_realigns:
	    if aligner == 'gmap':
		return_code = gmap.run(out_file, realign_bam_file, genome, index_dir, num_procs, multi=True)
	    elif aligner == 'bwa_mem':
		return_code = bwa_mem.run(out_file, realign_bam_file, genome, index_dir, num_procs)
	
	return realign_bam_file

	