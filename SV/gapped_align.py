from optparse import OptionParser
import re
from sets import Set
import pysam
from shared.alignment import Alignment, reverse_complement
from variant import Adjacency
from shared import gmap, bwa_mem

def find_single_unique(alns, aligner, bam, debug=False):
    return {
        'gmap': gmap.find_single_unique,
        'bwa_mem': bwa_mem.find_single_unique,
    }[aligner](alns, bam, debug=debug)
    
def find_adjs(align, contig_seq, is_transcriptome, ins_as_ins=False):
    adjs = []
    
    if not align.is_valid() or\
       not align.cigarstring or\
       len(align.blocks) != len(align.query_blocks) or\
       not re.search('[ID]\d+', align.cigarstring):
	return adjs
                
    # identify all (target) gaps in cigar: N=intron, D=deletion
    count = 0
    target_gaps = {}
    for gap in re.findall('[ND]\d+', align.cigarstring):
        target_gaps[count] = gap
        count += 1
        
    gap_count = 0
    for i in range(0, len(align.blocks) - 1):
        target_gap = align.blocks[i + 1][0] - align.blocks[i][1] - 1
        if align.strand == '+':
            query_gap = align.query_blocks[i + 1][0] - align.query_blocks[i][1] - 1
        else:
            query_gap = align.query_blocks[i][1] - align.query_blocks[i + 1][0] - 1
            
        # deletion or indel
        if target_gap > 0:
            if target_gaps.has_key(gap_count):
		breaks = (align.blocks[i][1], align.blocks[i + 1][0])
		contig_breaks = (align.query_blocks[i][1], align.query_blocks[i + 1][0])
				
		if target_gaps[gap_count][0] == 'N' and is_transcriptome:
		    continue
		
		if query_gap == 0:
		    rearrangement = 'del'
		else:
		    rearrangement = 'indel'

		probe_seq, break_pos = Adjacency.extract_probe(contig_seq, contig_breaks)
		#probe_seq = Adjacency.extract_probe(contig_seq, contig_breaks, kmer_size=51, min_buffer=4)
				
		adj = Adjacency((align.target, align.target), breaks, rearrangement, 
		                contig=align.query, contig_sizes=len(contig_seq), 
		                contig_breaks=contig_breaks, probes=probe_seq, orients=('L', 'R'),
		                aligns=[align],
		                align_types='gapped')
		adjs.append(adj)
	    
		
            gap_count += 1
            
        # insertion
        elif query_gap > 0:
	    event = 'ins'
	    breaks = (align.blocks[i][1], align.blocks[i][1])
	    contig_breaks = (align.query_blocks[i][1], align.query_blocks[i + 1][0])
	    
            if align.strand == '+':
                novel_seq = contig_seq[align.query_blocks[i][1] : align.query_blocks[i][1] + query_gap]
		novel_seq_ref = novel_seq
            else:
                novel_seq = contig_seq[align.query_blocks[i + 1][0] : align.query_blocks[i + 1][0] + query_gap]
		novel_seq_ref = reverse_complement(novel_seq)
		
	    length_novel_seq = len(novel_seq)
			    
	    duplicated = is_duplicated(novel_seq, sorted(contig_breaks, key=int), contig_seq)
	    if duplicated[0] > 0 or duplicated[1] > 0:
		if not ins_as_ins:
		    event = 'dup'
		
		if duplicated[0] > 0:
		    if align.strand == '+':
			breaks = (breaks[0] - length_novel_seq, breaks[1] + 1)
		    else:
			breaks = (breaks[0], breaks[0] + length_novel_seq + 1)
		else:
		    if align.strand == '+':
			breaks = (breaks[0], breaks[0] + length_novel_seq + 1)
		    else:
			breaks = (breaks[0] - length_novel_seq, breaks[1] + 1)
		    
		contig_breaks = (contig_breaks[0] - length_novel_seq * duplicated[0],
	                         contig_breaks[1] + length_novel_seq * duplicated[1])
		    			
	    probe_seq, break_pos = Adjacency.extract_probe(contig_seq, contig_breaks)
	    
            adj = Adjacency((align.target, align.target), breaks, event, contig=align.query, 
	                    contig_breaks=contig_breaks, contig_sizes=len(contig_seq), probes=probe_seq, 
	                    novel_seq=novel_seq_ref, orients=('L', 'R'), aligns=[align], align_types='gapped')
	    adjs.append(adj)
            
    return adjs

def is_duplicated(novel_seq, contig_breaks, contig_seq, min_len=3):    
    duplicated = [0, 0]
    length_novel_seq = len(novel_seq)
    if length_novel_seq > min_len:
	num_repeats = 0	
	# start = index position
	start = contig_breaks[0] - length_novel_seq
	while start >= 0:
	    before_seq = contig_seq[start : start + length_novel_seq]
	    if novel_seq.upper() == before_seq.upper():
		num_repeats += 1
		start -= length_novel_seq
	    else:
		break
	duplicated[0] = num_repeats
	
	num_repeats = 0	
	start = contig_breaks[1] - 1
	while start + length_novel_seq < len(contig_seq):
	    after_seq = contig_seq[start : start + length_novel_seq]
	    if novel_seq.upper() == after_seq.upper():
		num_repeats += 1
		start += length_novel_seq
	    else:
		break
	duplicated[1] = num_repeats
		
    return duplicated
	    
def screen_probe_alns(adj_aligns, probe_alns, align_type, min_pc_mapped=1.0):
    for aln in probe_alns:
	matched_len = 0
	query_len = sum([a[1] for a in aln.cigar if a[0] in (0, 1, 4, 5)])
	matched_len = sum([a[1] for a in aln.cigar if a[0] == 0])
	# allows gap inside probe's alignment for checking end-to-end mapping
	if re.match('^\d+M.+\d+M$', aln.cigarstring):
	    matched_len = query_len
		
	#if align_type == 'split' and (re.match('\d+M$', aln.cigarstring) or re.match('\d+M.+\d+M$', aln.cigarstring)):
	if align_type == 'split' and float(matched_len) / float(query_len) >= min_pc_mapped:
	    return False
	
	# if it's a single alignment and the probe can map perfectly (no clips, no insertion/deletion) to a location -> out
	if align_type == 'gapped' and aln.rlen == aln.alen and not re.search('[DIN]', aln.cigarstring):
	    return False
	
    return True
    