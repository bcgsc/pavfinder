from sets import Set
from SV.split_align import call_event
from SV.variant import Adjacency

class FusionFinder:
    @classmethod
    def find_chimera(cls, chimera_block_matches, transcripts, aligns):
	"""Identify gene fusion between split alignments
	
	Args:
	    chimera_block_matches: (list) dictionaries where 
	                                      key=transcript name, 
	                                      value=[match1, match2, ...] where
						    match1 = matches of each alignment block
							     i.e.
							     [(exon_id, '=='), (exon_id, '==')] 
	    transcripts: (dict) key=transcript_name value=Transcript object
	    aligns: (list) Alignments involved in chimera
	Returns:
	    Adjacency with genes, transcripts, exons annotated
	"""
	assert len(chimera_block_matches) == len(aligns), 'number of block matches(%d) != number of aligns(%d)' % \
	       (len(chimera_block_matches), len(aligns))
	for (matches1_by_txt, matches2_by_txt) in zip(chimera_block_matches, chimera_block_matches[1:]):	    
	    genes1 = Set([transcripts[txt].gene for txt in matches1_by_txt.keys()])
	    genes2 = Set([transcripts[txt].gene for txt in matches2_by_txt.keys()])
	    	    	    
	    junc_matches1 = {}
	    num_blocks = len(aligns[0].blocks)
	    for transcript in chimera_block_matches[0].keys():
		junc_matches1[transcript] = chimera_block_matches[0][transcript][num_blocks - 1]
	
	    junc_matches2 = {}
	    for transcript in chimera_block_matches[1].keys():
		junc_matches2[transcript] = chimera_block_matches[1][transcript][0]
    
	    # create adjacency first to establish L/R orientations, which is necessary to pick best transcripts
	    fusion = call_event(aligns[0], aligns[1], no_sort=True)
	    junc1, junc2 = cls.identify_fusion(junc_matches1, junc_matches2, transcripts, fusion.orients)
	    if junc1 and junc2:
		cls.annotate_fusion(fusion, junc1, junc2, transcripts)
		return fusion
	    
	return None

    @classmethod
    def find_read_through(cls, matches_by_transcript, transcripts, align):
	"""Identify read-through fusion
	
	Assume single alignment is spanning 2 genes
	
	Args:
	    matches_by_transcript: (list) dictionaries where 
	                                      key=transcript name, 
	                                      value=[match1, match2, ...] where
						    match1 = matches of each alignment block
							     i.e.
							     [(exon_id, '=='), (exon_id, '==')] 
							     or None if no_match
	    transcripts: (dict) key=transcript_name value=Transcript object
	    align: (Alignment) alignment spanning >1 gene
	Returns:
	    Adjacency with genes, transcripts, exons annotated
	"""
	def create_fusion(junc1, junc2, pos, contig_breaks):
	    """Creates Adjacency object capturing read-through
	    
	    Args:
	        junc1: (tuple) transcript_id, 'match' list [(exon_id, '==')]
		junc2: (tuple) transcript_id, 'match' list [(exon_id, '==')]
		pos: (tuple) genome coordinates of fusion
		contig_breaks: (tuple) contig coordinates of breaks
	    Returns:
		Adjacency with genes, transcripts, exons annotated
	    """
	    fusion = Adjacency((align.target, align.target),
	                       pos,
	                       '-',
	                       orients=('L', 'R'),
	                       contig_breaks = contig_breaks
	                       )     
	    cls.annotate_fusion(fusion, junc1, junc2, transcripts)
	    return fusion
	
	num_blocks = len(align.blocks)
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
	for k in range(len(junctions)):
	    i, j = junctions[k]
	    if len(genes_in_block[i]) == 1 and len(genes_in_block[j]) == 1:
		if list(genes_in_block[i])[0] != list(genes_in_block[j])[0]:
		    junc1, junc2 = cls.identify_fusion(matches_by_block[i], matches_by_block[j], transcripts, ('L', 'R'))
		    pos = (align.blocks[i][1], align.blocks[j][0])
		    contig_breaks = (align.query_blocks[i][1], align.query_blocks[j][0])
		    fusion = create_fusion(junc1, junc2, pos, contig_breaks)
		    return fusion
	
	    elif len(genes_in_block[i]) > 1:
		print 'fusion block', i, j, genes_in_block
		
	    elif len(genes_in_block[j]) > 1:
		print 'fusion block', i, j
		
	for i in genes_in_block.keys():
	    if len(genes_in_block[i]) == 2:
		cls.identify_fusion_unknown_break(matches_by_block[i], transcripts)
		
	return None
		
    @classmethod
    def identify_fusion_unknown_break(cls, matches, transcripts):
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
    
    @classmethod
    def identify_fusion(cls, matches1, matches2, transcripts, orients):
	"""Given 2 block matches pick the 2 transcripts"""
	def pick_best(matches, orient, exon_bound_score=1000):
	    scores = {}
	    for transcript in matches:
		if orient == 'L':
		    junction_block = -1
		    junction_edge, distant_edge = -1, 0
		else:
		    junction_block = 0
		    junction_edge, distant_edge = 0, -1
		score = 0
		if matches[transcript] is not None:
		    if matches[transcript][junction_block][1][junction_edge] == '=':
			score += exon_bound_score
			
		    # align block within exon gets more points
		    if matches[transcript][junction_block][1][distant_edge] == '=':
			score += 15
		    elif matches[transcript][junction_block][1][distant_edge] == '>':
			score += 10
		    elif matches[transcript][junction_block][1][distant_edge] == '<':
			score += 5
		scores[transcript] = score
		
	    best_score = max(scores.values())
	    best_txt = sorted([t for t in matches.keys() if scores[t] == best_score], 
		               key=lambda t: transcripts[t].length(), reverse=True)[0]
	    return best_txt, best_score > exon_bound_score
	    	
	best_txt1, exon_bound1 = pick_best(matches1, orients[0])
	best_txt2, exon_bound2 = pick_best(matches2, orients[1])
	
	return (best_txt1, matches1[best_txt1], exon_bound1), (best_txt2, matches2[best_txt2], exon_bound2)

    @classmethod
    def is_ptd(cls, fusion):
	"""Define PTD"""
	if fusion.genes[0] == fusion.genes[1]:
	    fusion.rna_event = 'PTD'
	    
    @classmethod
    def annotate_fusion(cls, fusion, junc1, junc2, transcripts):
	fusion.rna_event = 'fusion'
	fusion.genes = (transcripts[junc1[0]].gene, transcripts[junc2[0]].gene)
	fusion.transcripts = (junc1[0], junc2[0])
	#if junc1[1] is None or junc2[1] is None:
	    #print 'check', junc1[0], junc2[0], transcripts[junc1[0]], transcripts[junc2[0]], junc1, junc2
	    #return None
	fusion.exons = transcripts[junc1[0]].exon_num(junc1[1][0][0]), transcripts[junc2[0]].exon_num(junc2[1][0][0])
	fusion.exon_bound = junc1[2], junc2[2]
	fusion.in_frame = True if (fusion.exon_bound[0] and fusion.exon_bound[1]) else False
	
	sense_fusion = cls.is_sense(fusion, transcripts[junc1[0]], transcripts[junc2[0]], fusion.orients[0], fusion.orients[1])
	if not sense_fusion:
	    fusion.is_sense = False
	    fusion.gene5, fusion.gene3 = '-', '-'
	else:
	    fusion.is_sense = True
	    fusion.gene5, fusion.gene3 = sense_fusion[0].gene, sense_fusion[1].gene
	
	cls.is_ptd(fusion)
	    
    @classmethod
    def is_sense(cls, fusion, transcript1, transcript2, orient1, orient2):
	t5, t3 = None, None
	if fusion.exon_bound[0] and fusion.exon_bound[1]:
	    if orient1 == 'L' and transcript1.strand == '+' and orient2 == 'R' and transcript2.strand == '+':
		t5, t3 = transcript1, transcript2
	    elif orient1 == 'L' and transcript1.strand == '-' and orient2 == 'R' and transcript2.strand == '-':
		t3, t5 = transcript1, transcript2
	    elif orient1 == 'R' and transcript1.strand == '+' and orient2 == 'L' and transcript2.strand == '+':
		t3, t5 = transcript1, transcript2
	    elif orient1 == 'R' and transcript1.strand == '-' and orient2 == 'L' and transcript2.strand == '-':
		t5, t3 = transcript1, transcript2
		
	print 'abc', transcript1.gene, transcript2.gene, orient1, orient2, fusion.exon_bound, t5, t3
	if t5 is not None and t3 is not None:
	    return (t5, t3)
	else:
	    return False
