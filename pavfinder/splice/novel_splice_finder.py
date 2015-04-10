from sets import Set
from pavfinder.shared.alignment import reverse_complement

def find_novel_junctions(block_matches, align, transcripts, ref_fasta, accessory_known_features):
    """Find novel junctions within a single gene/transcript
    
    Args:
	block_matches: (list) dictionaries where 
					  key=transcript name, 
					  value=[match1, match2, ...] where
						match1 = matches of each alignment block
							 i.e.
							 [(exon_id, '=='), (exon_id, '==')] 
	align: (Alignment) alignment object
	transcripts: (dictionary) key=transcript_id value=Transcript object
	ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
    Returns:
	List of event (dictionary storing type of event, exon indices, genomic coordinates and size)
    """
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
	    
	    if is_junction_annotated(matches[i][-1], matches[j][0]):
		annotated.add((i, j))
		continue
	    
    known_exons = known_jns = None
    if accessory_known_features is not None:
	known_exons = accessory_known_features['exon']
	known_jns = accessory_known_features['junction']

    all_events = []
    for transcript in block_matches.keys():
	matches = block_matches[transcript]
	for i in range(len(matches) - 1):
	    j = i + 1
		    
	    if matches[i] is None and matches[j] is not None:
		# special case where the first 2 blocks is the utr and there's an insertion separating the 2 blocks
		if i == 0:
		    events = classify_novel_junction(matches[i], 
		                                     matches[j][0], 
		                                     align.target, 
		                                     align.blocks[i:j+1], 
		                                     transcripts[transcript],
		                                     ref_fasta,
		                                     known_jns = known_jns,
		                                     known_exons = known_exons,
		                                     )
		    for e in events:
			e['blocks'] = (i, j)
			e['transcript'] = transcript
			e['contig_breaks'] = (align.query_blocks[i][1], align.query_blocks[j][0])
		    all_events.extend(events)
		continue
	    
	    # for retained intron, not a 'junction'
	    if matches[i] is not None and len(matches[i]) > 1:
		events = classify_novel_junction(matches[i], 
		                                 None,
		                                 align.target, 
		                                 align.blocks[i], 
		                                 transcripts[transcript],
		                                 ref_fasta,
		                                 known_jns = known_jns,
		                                 known_exons = known_exons,
		                                 )
		if events:
		    for e in events:
			e['blocks'] = (i, j)
			e['transcript'] = transcript
			e['contig_breaks'] = (align.query_blocks[i][0], align.query_blocks[i][1])
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
		events = classify_novel_junction(matches[i][-1], 
		                                 matches[j][0], 
		                                 align.target, 
		                                 align.blocks[i:j+1], 
		                                 transcripts[transcript],
		                                 ref_fasta,
		                                 known_jns = known_jns,
		                                 known_exons = known_exons,
		                                 )
		if events:
			for e in events:
			    e['blocks'] = range(i + 1, j)
			    e['contig_breaks'] = (align.query_blocks[i][1], align.query_blocks[j][0])
			    e['transcript'] = transcript
			all_events.extend(events)
	    
    if all_events:
	# group events
	grouped = {}
	for event in all_events:
	    key = str(event['blocks']) + event['event']
	    try:
		grouped[key].append(event)
	    except:
		grouped[key] = [event]
	    
	uniq_events = []
	for events in grouped.values():
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
	                        'contig_breaks': [contig_breaks],
	                        'size': events[0]['size'],
	                        }
	                       )
		
	return uniq_events
    
def classify_novel_junction(match1, match2, chrom, blocks, transcript, ref_fasta, min_intron_size=20, known_jns=None, known_exons=None):
    """Classify given junction into different splicing events or indel
    
    Args:
	match1: (tuple or list) single tuple: exon_index, 2-char match result e.g. '==', '>=', etc
				list of 2 tuples for retained_intron (special case)
	match2: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
	chrom: (str) chromosome, for checking splice motif
	blocks: (list) list of list (block coordinates)
	transcript: (Transcript) object of transcript, for getting exon coordinates
	ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
    Returns:
	List of event (dictionary storing type of event, exon indices, genomic coordinates and size)
    """
    events = []
    # set default values for 'pos'
    if type(blocks[0]) is int:
	pos = (blocks[0], blocks[1])
	if known_exons and (chrom, pos[0], pos[1]) in known_exons:
	    return events
    else:
	pos = (blocks[0][1], blocks[1][0])
	if known_jns and (chrom, pos[0], pos[1]) in known_jns:
	    return events

    if match2 is None:
	if len(match1) == 2:
	    exons = [m[0] for m in match1]
	    if match1[0][1] == '=>' and\
	       match1[-1][1] == '<=' and\
	       len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(match1) - 1:
		size = transcript.exons[exons[1]][0] - transcript.exons[exons[0]][1] - 1 
		events.append({'event': 'retained_intron', 'exons': exons, 'pos':pos, 'size':size})
	      
    # genomic deletion or novel_intron
    elif match1 is None and match2 is not None:
	if match2[0] == 0 and match2[1] == '<=':
	    if transcript.strand == '+':
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
		
	    gap_size = blocks[1][0] - blocks[0][1] - 1
	    pos = (blocks[0][1], blocks[1][0])
	    event = None
	    if gap_size > 0:
		if gap_size > min_intron_size and\
	           check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		    event = 'novel_intron'
		else:
		    event = 'del'
	    if event is not None:
		events.append({'event': event, 'exons': [match2[0]], 'pos':pos, 'size':gap_size})
	    
    else:
	if match2[0] > match1[0] + 1 and\
           '=' in match1[1] and\
           '=' in match2[1]:
	    size = 0
	    for e in range(match1[0] + 1, match2[0]):
		exon = transcript.exons[e]
		size += exon[1] - exon[0] + 1
	    events.append({'event': 'skipped_exon', 'exons': range(match1[0] + 1, match2[0]), 'pos':pos, 'size':size})
	   
	if match1[0] == match2[0] and\
           match1[1][1] == '<' and match2[1][0] == '>':		
	    if transcript.strand == '+':
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
	    
	    gap_size = blocks[1][0] - blocks[0][1] - 1
	    pos = (blocks[0][1], blocks[1][0])
	    event = None
	    if gap_size > 0:
		if gap_size > min_intron_size and\
	           check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		    event = 'novel_intron'
		else:
		    event = 'del'
	    elif gap_size == 0:
		event = 'ins'
	    
	    if event is not None:
		if event != 'ins':
		    events.append({'event': event, 'exons': [match1[0]], 'pos':pos, 'size':gap_size})
		else:
		    events.append({'event': event, 'exons': [match1[0]], 'pos':pos})
	    
	# novel donor and acceptor
	if match1[1][1] == '=' and match2[1][0] != '=':
	    size = abs(blocks[1][0] - transcript.exons[match2[0]][0])
	    if transcript.strand == '+':
		event = 'novel_acceptor'
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		event = 'novel_donor'
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
	    # check splice motif
	    if check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		events.append({'event': event, 'exons': [match1[0], match2[0]], 'pos':pos, 'size':size})
		
	if match1[1][1] != '=' and match2[1][0] == '=':
	    size = abs(blocks[0][1] - transcript.exons[match1[0]][1])
	    if transcript.strand == '+':
		event = 'novel_donor'
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		event = 'novel_acceptor'
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
	    # check splice motif
	    if check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		events.append({'event': event, 'exons': [match1[0], match2[0]], 'pos':pos, 'size':size})
		    
	if match2[0] == match1[0] + 1 and\
           match1[1][1] == '=' and\
           match2[1][0] == '=': 
	    pos = (blocks[1][0], blocks[-2][1])
	    size = blocks[-2][1] - blocks[1][0] + 1
	    if transcript.strand == '+':
		donor_start = pos[1] + 1
		acceptor_start = pos[0] - 2
	    else:
		donor_start = pos[0] - 2
		acceptor_start = pos[1] + 1
	    # check splice motif
	    if check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		events.append({'event': 'novel_exon', 'exons': [], 'pos':pos, 'size':size})
	    
    # set size to None for event that doesn't have size i.e. 'ins'
    for event in events:
	if not event.has_key('size'):
	    event['size'] = None
    
    return events

def check_splice_motif(chrom, donor_start, acceptor_start, strand, ref_fasta):
    """Check if the 4-base splice motif of a novel junction is canonical (gtag)
	
    Right now only considers 'gtag' as canonical
    
    Args:
	chrom: (str) chromosome
	donor_start: (int) genomic position of first(smallest) base of donor site (1-based)
	acceptor_start: (int) genomic position of first(smallest) base of acceptor site (1-based)
	strand: (str) transcript strand '+' or '-'
	ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence
    Returns:
	True if it's canonical, False otherwise
    """
    canonical_motifs = Set()
    canonical_motifs.add('gtag')
    
    donor_seq = ref_fasta.fetch(chrom, donor_start - 1, donor_start - 1 + 2)
    acceptor_seq = ref_fasta.fetch(chrom, acceptor_start - 1, acceptor_start - 1 + 2)
    
    if strand == '+':
	motif = donor_seq + acceptor_seq
    else:
	motif = reverse_complement(acceptor_seq + donor_seq)
    
    if motif.lower() in canonical_motifs:
	return True
    else:
	return False
    
def is_junction_annotated(match1, match2):
    """Checks if junction is in gene model
    
    Args:
	match1: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
	match2: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
    Returns:
	True or False
    """
    if match2[0] == match1[0] + 1 and\
       match1[1][1] == '=' and\
       match2[1][0] == '=':
	return True
    
    return False

def annotate_ref_junctions(events, junction_depths, transcripts):
    """Annotate 5' and 3' gene reference junction depth/coverage

    Arguments:
        events: (list) of events (Event)
        juncton_depths: (dict) {chrom[(start, end)] = depth}
        transcripts: (dict) transcript ID to object mapping (for extracting exon coordinates)
    """
    for event in events:
	txt = transcripts[event.transcripts[0]]
	jn5 = None
	jn3 = None
	if event.rna_event == 'novel_donor':
	    jn5 = txt.exon(event.exons[0])[1], txt.exon(event.exons[1])[0]

	elif event.rna_event == 'novel_acceptor':
	    jn5 = txt.exon(event.exons[0])[1], txt.exon(event.exons[1])[0]

	elif event.rna_event == 'skipped_exon':
	    exon5 = min(event.exons) - 1
	    exon3 = max(event.exons) + 1
	    exon5_coord = txt.exon(exon5)
	    exon3_coord = txt.exon(exon3)
	    if exon5_coord is not None:
		exon_coord = txt.exon(min(event.exons))
		if txt.strand == '+':
		    jn5 = exon5_coord[1], exon_coord[0]
		else:
		    jn5 = exon_coord[1], exon5_coord[0]

	    if exon3_coord is not None:
		exon_coord = txt.exon(max(event.exons))
		if txt.strand == '+':
		    jn3 = exon_coord[1], exon3_coord[0]
		else:
		    jn3 = exon3_coord[1], exon_coord[0]

	elif event.rna_event == 'retained_intron':
	    jn5 = txt.exon(event.exons[0])[1], txt.exon(event.exons[1])[0]

	elif event.rna_event == 'novel_exon':
	    for i in range(len(txt.exons) - 1):
		intron = txt.exons[i][1] + 1, txt.exons[i + 1][0] - 1
		if event.breaks[0] > intron[0] and event.breaks[1] < intron[1]:
		    jn5 = txt.exons[i][1], txt.exons[i + 1][0]
		    break

	elif event.rna_event == 'novel_intron':
	    exon = event.exons[0]
	    exon_coord = txt.exon(exon)
	    exon5 = exon - 1
	    exon3 = exon + 1
	    exon5_coord = txt.exon(exon5)
	    exon3_coord = txt.exon(exon3)
	    if exon5_coord is not None:
		if txt.strand == '+':
		    jn5 = exon5_coord[1], exon_coord[0]
		else:
		    jn5 = exon_coord[1], exon5_coord[0]

	    if exon3_coord is not None:
		if txt.strand == '+':
		    jn3 = exon_coord[1], exon3_coord[0]
		else:
		    jn3 = exon3_coord[1], exon_coord[0]

	if jn5 is not None or jn3 is not None:
	    if junction_depths.has_key(event.chroms[0]):
		if jn5 is not None and junction_depths[event.chroms[0]].has_key(jn5):
		    event.ref5_depth = junction_depths[event.chroms[0]][jn5]
		    event.ref5_coord = '%s:%s-%s' % (event.chroms[0], jn5[0], jn5[1])
		if jn3 is not None and junction_depths[event.chroms[0]].has_key(jn3):
		    event.ref3_depth = junction_depths[event.chroms[0]][jn3]
		    event.ref3_coord = '%s:%s-%s' % (event.chroms[0], jn3[0], jn3[1])
