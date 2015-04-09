from sets import Set
from itertools import groupby
import sys
import re
import pysam
from pavfinder.SV.split_align import call_event
from pavfinder.shared.adjacency import Adjacency
from pavfinder.SV import gapped_align, split_align
from pavfinder.splice.event import Event

def find_chimera(chimera_block_matches, transcripts, aligns, contig_seq, exon_bound_only=False, coding_only=False, sense_only=False):
    """Identify gene fusion between split alignments
	
    Args:
        chimera_block_matches: (list) dictionaries where
                                          key=transcript name,
                                          value=[match1, match2, ...] where
                                                match1 = matches of each alignment block
                                                         i.e.
                                                         [(exon_id, '=='), (exon_id, '==')]
        transcripts: (dict) key=transcript_name value=Transcript object
        aligns: (list) ORDERED alignments involved in chimera
    Returns:
        Adjacency with genes, transcripts, exons annotated
    """
    assert len(chimera_block_matches) == len(aligns), 'number of block matches(%d) != number of aligns(%d)' % \
           (len(chimera_block_matches), len(aligns))
    for (matches1_by_txt, matches2_by_txt) in zip(chimera_block_matches, chimera_block_matches[1:]):
	# create adjacency first to establish L/R orientations
	# necessary for choosing block_matches(start or end) and picking best transcripts
	fusion = call_event(aligns[0], aligns[1], no_sort=True, contig_seq=contig_seq, probe_side_len=50)
	fusion.__class__ = Event
	fusion.read_depth = 0

	genes1 = Set([transcripts[txt].gene for txt in matches1_by_txt.keys()])
	genes2 = Set([transcripts[txt].gene for txt in matches2_by_txt.keys()])

	junc_matches1 = {}
	for transcript in chimera_block_matches[0].keys():
	    if fusion.orients[0] == 'L':
		junc_matches1[transcript] = chimera_block_matches[0][transcript][-1]
	    else:
		junc_matches1[transcript] = chimera_block_matches[0][transcript][0]

	junc_matches2 = {}
	for transcript in chimera_block_matches[1].keys():
	    if fusion.orients[1] == 'L':
		junc_matches2[transcript] = chimera_block_matches[1][transcript][-1]
	    else:
		junc_matches2[transcript] = chimera_block_matches[1][transcript][0]

	junc1, junc2 = identify_fusion(junc_matches1, junc_matches2, transcripts, fusion.orients)
	if junc1 and junc2 and junc1[1] is not None and junc2[1] is not None:
	    annotate_fusion(fusion, junc1, junc2, transcripts)
	    if exon_bound_only and not (fusion.exon_bound[0] and fusion.exon_bound[1]):
		return None

	    if coding_only and not (transcripts[junc1[0]].coding and transcripts[junc2[0]].coding):
		return None

	    if sense_only and not fusion.is_sense:
		return None

	    return fusion
	
    return None

def find_read_through(matches_by_transcript, transcripts, align, exon_bound_only=False, coding_only=False, sense_only=False):
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
	fusion = Event((align.target, align.target),
                       pos,
                       '-',
                       orients=('L', 'R'),
                       contig_breaks = contig_breaks
                       )

	if exon_bound_only and not (fusion.exon_bound[0] and fusion.exon_bound[1]):
	    return None
    
	if coding_only and not (transcripts[junc1[0]].coding and transcripts[junc2[0]].coding):
	    return None
	
	if sense_only and not fusion.is_sense:
	    return None
	
	annotate_fusion(fusion, junc1, junc2, transcripts)
	fusion.rna_event = 'read_through'
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
		junc1, junc2 = identify_fusion(matches_by_block[i], matches_by_block[j], transcripts, ('L', 'R'))
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
	    identify_fusion_unknown_break(matches_by_block[i], transcripts)

    return None

def identify_fusion_unknown_break(matches, transcripts):
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

def identify_fusion(matches1, matches2, transcripts, orients):
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

def is_ptd(fusion):
    """Define PTD"""
    if fusion.genes[0] == fusion.genes[1] and fusion.exon_bound[0] and fusion.exon_bound[1]:
	fusion.rna_event = 'PTD'
	
def annotate_fusion(fusion, junc1, junc2, transcripts):
    fusion.rna_event = 'fusion'
    fusion.genes = (transcripts[junc1[0]].gene, transcripts[junc2[0]].gene)
    fusion.transcripts = (junc1[0], junc2[0])
    fusion.exons = transcripts[junc1[0]].exon_num(junc1[1][0][0]), transcripts[junc2[0]].exon_num(junc2[1][0][0])
    fusion.exon_bound = junc1[2], junc2[2]
    fusion.in_frame = True if (fusion.exon_bound[0] and fusion.exon_bound[1]) else False
    fusion.ref5_depth = fusion.ref3_depth = fusion.ref5_coord = fusion.ref3_coord = None
	    
    sense_fusion = is_sense(fusion, transcripts[junc1[0]], transcripts[junc2[0]], fusion.orients[0], fusion.orients[1])
    if not sense_fusion:
	fusion.is_sense = False
	fusion.gene5, fusion.gene3 = '-', '-'
    else:
	fusion.is_sense = True
	fusion.gene5, fusion.gene3 = sense_fusion[0].gene, sense_fusion[1].gene

    is_ptd(fusion)

    # chimera within same gene but not PTD, change event label to rearrangement instead of 'fusion'
    if fusion.genes[0] == fusion.genes[1] and fusion.rna_event != 'PTD':
	fusion.rna_event = fusion.rearrangement

def is_sense(fusion, transcript1, transcript2, orient1, orient2):
    t5, t3 = None, None
    if (orient1 == 'L' and transcript1.strand == '+' and orient2 == 'R' and transcript2.strand == '+') or\
       (orient1 == 'R' and transcript1.strand == '-' and orient2 == 'L' and transcript2.strand == '-') or\
       (orient1 == 'L' and transcript1.strand == '+' and orient2 == 'L' and transcript2.strand == '-') or\
       (orient1 == 'R' and transcript1.strand == '-' and orient2 == 'R' and transcript2.strand == '+'):
	t5, t3 = transcript1, transcript2
    elif (orient1 == 'L' and transcript1.strand == '-' and orient2 == 'R' and transcript2.strand == '-') or\
         (orient1 == 'R' and transcript1.strand == '+' and orient2 == 'L' and transcript2.strand == '+') or\
         (orient1 == 'L' and transcript1.strand == '-' and orient2 == 'L' and transcript2.strand == '+') or\
         (orient1 == 'R' and transcript1.strand == '+' and orient2 == 'R' and transcript2.strand == '-'):
	t3, t5 = transcript1, transcript2
	
	    
    if t5 is not None and t3 is not None:
	return (t5, t3)
    else:
	return False

def screen_realigns(fusions, outdir, align_info, contigs_fasta=None, name_sep='-', debug=False):
    """Screens fusions by religning probes against reference genome

    Will filter out events whose probe can align to single location

    Args:
	fusions: (list) Adjacencies
	outdir: (str) absolute path of output directory, for storing re-alignment results
	aligner: (str) aligner name (gmap)
	align_info: (dict) 'genome', 'index_dir', 'num_procs'
	name_sep: (str) single character to separate event info in generating query name
	debug: (boolean) output debug info e.g. reason for screening out event
    """
    if align_info is not None:
	realign_bam_file = Adjacency.realign(fusions, outdir, align_info['aligner'],
                                             probe=True, subseq=True,
                                             contigs_fasta=contigs_fasta,
                                             name_sep=name_sep,
                                             num_procs=align_info['num_procs'],
                                             genome=align_info['genome'],
                                             index_dir=align_info['index_dir'])
	
    try:
	bam = pysam.Samfile(realign_bam_file, 'rb')
    except:
	sys.exit('Error parsing realignment BAM:%s' % realign_bam_file)
	
    # creates mapping from query to fusion
    query_to_fusion = {}
    for i in range(len(fusions)):
	query = fusions[i].contigs[0] + name_sep + fusions[i].key()
	query_to_fusion[query] = i
	    
    # captures contig names whose probe aligns to single position
    # will remove any events coming from these contigs because the split-alignment is not reliable
    failed_contigs = Set()
    for key, group in groupby(bam.fetch(until_eof=True), lambda x: name_sep.join(x.qname.split(name_sep)[:8])):
	alns = list(group)
	fusion_idx = query_to_fusion[key]
	fusion = fusions[fusion_idx]
	fusion_aligns = fusion.aligns[0]
		
	probe_alns = [aln for aln in alns if not aln.qname[-1].isdigit()]
	if not gapped_align.screen_probe_alns(fusion_aligns, probe_alns, fusion.align_types[0], min_pc_mapped=0.9):
	    if debug:
		sys.stdout.write('probe align completely to one location: %s - contigs(s) filtered out\n' % key)
	    for contig in fusion.contigs:
		failed_contigs.add(contig)
	    continue
	    
	subseq_alns = [aln for aln in alns if aln.qname[-1].isdigit()]
	specific1, specific2 = split_align.screen_subseq_alns(fusion_aligns, subseq_alns, bam, name_sep, debug=debug)
	if not specific1 or not specific2:
	    for contig in fusion.contigs:
		failed_contigs.add(contig)
	    continue

    return failed_contigs

def annotate_ref_junctions(fusions, junction_depths, transcripts):
    """Annotate 5' and 3' gene reference junction depth/coverage

    Arguments:
        fusions: (list) of fusions (Event)
        juncton_depths: (dict) {chrom[(start, end)] = depth}
        transcripts: (dict) transcript ID to object mapping (for extracting exon coordinates)
    """
    def update_depth(fusion, depth, gene, coord):
	"""Report the highest-expressed reference junction"""
	if gene == fusion.gene5:
	    if fusion.ref5_depth is None or depth > fusion.ref5_depth:
		fusion.ref5_depth = depth
		fusion.ref5_coord = coord
	elif gene == fusion.gene3:
	    if fusion.ref3_depth is None or depth > fusion.ref3_depth:
		fusion.ref3_depth = depth
		fusion.ref3_coord = coord

    for fusion in fusions:
	txts = (transcripts[fusion.transcripts[0]], transcripts[fusion.transcripts[1]])
	for i in (0,1):
	    # exon coordinates
	    jn_exon = txts[i].exon(fusion.exons[i])
	    if jn_exon is None:
		print 'cannot find junction', fusion.contigs, fusion.exons[i]
		continue
	    if fusion.orients[i] == 'L':
		jn_pos = jn_exon[1]
	    else:
		jn_pos = jn_exon[0]

	    for pos, depth in junction_depths[fusion.chroms[i]].iteritems():
		if fusion.orients[i] == 'L' and int(pos[0]) == int(jn_pos):
		    update_depth(fusion, depth, fusion.genes[i], '%s:%s-%s' % (fusion.chroms[i], pos[0], pos[1]))
		elif fusion.orients[i] == 'R' and int(pos[1]) == int(jn_pos):
		    update_depth(fusion, depth, fusion.genes[i], '%s:%s-%s' % (fusion.chroms[i], pos[0], pos[1]))
