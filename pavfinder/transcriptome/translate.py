from Bio.Seq import Seq
import re
import sys
from math import ceil

def get_orfs(sequence, frames=None, only_strand=None, min_size=1):
    def prot2nuc_coord(prot_coord, frame, strand, start=False):
        if prot_coord == 0:
            nuc_coord = 0
        else:
            if start:
                nuc_coord = prot_coord * 3
            else:
                nuc_coord = prot_coord * 3 + 2
        
        if strand == '+':
            nuc_coord += frame
        else:
            nuc_coord = len(sequence) - nuc_coord - 1 - frame

        return nuc_coord
    
    def get_seq(coords):
        if coords[1] < coords[0]:
            start, end = coords[1], coords[0]
        else:
            start, end = coords[0], coords[1]
        return sequence[start:end+1]
    
    all_orfs = []
    seq = Seq(sequence)
    seq_len = len(sequence)
    for strand, nuc in [('+', seq), ('-', seq.reverse_complement())]:
        if only_strand is not None and strand != only_strand:
            continue
        
        for frame in range(3):
            if frames is not None and not frame in frames:
                continue
            
            translation = str(nuc[frame:].translate())
            translation_len = len(translation)
            
            stops = [i.start() for i in re.finditer("\*", translation)]
            if stops:
                start = 0
                for stop in stops:
                    orf = translation[start:stop]
                    seq_start = prot2nuc_coord(start, frame, strand, start=True)
                    seq_end = prot2nuc_coord(stop - 1, frame, strand)
                    start = stop + 1
                    if len(orf) > min_size:
                        all_orfs.append((seq_start, seq_end, strand, frame, orf))
                if stop + 1 < len(translation):
                    start = stop + 1
                    stop = len(translation)
                    orf = translation[start:stop]
                    seq_start = prot2nuc_coord(start, frame, strand, start=True)
                    seq_end = prot2nuc_coord(stop - 1, frame, strand)
                    if len(orf) > min_size:
                        all_orfs.append((seq_start, seq_end, strand, frame, orf))
            else:
                start = 0
                stop = len(translation)
                orf = translation
                seq_start = prot2nuc_coord(start, frame, strand, start=True)
                seq_end = prot2nuc_coord(stop - 1, frame, strand)
                if len(orf) > min_size:
                    all_orfs.append((seq_start, seq_end, strand, frame, orf))

    orfs_sorted = sorted(all_orfs, key=lambda orf: len(orf[4]), reverse=True)
    return orfs_sorted

def check_frame(query_seq, seq_breaks, genome_breaks, genome_fasta, transcripts, event_type, size):
    in_frame = None
    if transcripts[0].is_coding() and transcripts[1].is_coding():
	if event_type in ('ins', 'del', 'repeat-expansion', 'repeat-reduction'):
	    if transcripts[0] == transcripts[1] and\
	       not transcripts[0].within_utr(genome_breaks[0]) and\
	       not transcripts[0].within_utr(genome_breaks[1]):
		if size % 3 == 0:
		    return True
		else:
		    return False
	else:
	    in_frame = is_inframe(transcripts[0],
		                  transcripts[1],
		                  seq_breaks,
		                  query_seq,
		                  genome_fasta)
    return in_frame

def is_inframe(txt5, txt3, query_breaks, query_seq, genome_fasta):
    """Checks if event is in frame
    
    Arguments:
        txt5: 5' Transcript
        txt3: 3' Transcript
        query_breaks: (tuple) query break coordinates (sorted)
        query_seq: (str) query sequence
        genome_fasta: (Fasta) genome Fasta for constructing transcript sequence
        
    Returns:
        tuple of (amino acid number of 5'transcript,
                  amino acid number of 3'transcript,
                  novel_amino acid in between)
    """
    max_aa_len_test = 5
    orf5 = get_orfs(txt5.get_sequence(genome_fasta, cds_only=True), frames=[0], only_strand='+')[0][-1]
    orf3 = get_orfs(txt3.get_sequence(genome_fasta, cds_only=True), frames=[0], only_strand='+')[0][-1]
    
    query_orfs = get_orfs(query_seq)
    for orf in query_orfs:
        orf_coords = sorted(orf[:2])
        
        if orf[2] == '-':
            orf_start, orf_end = len(query_seq) - orf[0], len(query_seq) - orf[1]
            break_start, break_end = len(query_seq) - query_breaks[1] + 1, len(query_seq) - query_breaks[0] + 1
        else:
            orf_start, orf_end = orf[0] + 1, orf[1] + 1
            break_start, break_end = query_breaks[0], query_breaks[1]
            
        if orf_start < break_start and orf_end > break_end:
            ctg_breaks_aa = int(ceil((break_start - orf_start) / 3.0)), int(ceil((break_end - orf_start) / 3.0)) + 1
            if ctg_breaks_aa[1] - ctg_breaks_aa[0] > 1:
                novel_aa = orf[-1][ctg_breaks_aa[0]:ctg_breaks_aa[1] - 1]
            else:
                novel_aa = None
                        
            up_orf = orf[-1][:ctg_breaks_aa[0]]
            down_orf = orf[-1][ctg_breaks_aa[1] - 1:]
            
            if len(up_orf) < max_aa_len_test or len(down_orf) < max_aa_len_test:
                continue
            
            up_orf_test = up_orf[-1 * min(max_aa_len_test, len(up_orf)):]
            down_orf_test = down_orf[:min(max_aa_len_test, len(down_orf))]
	    down_orf_test1 = down_orf[1:][:min(max_aa_len_test, len(down_orf))]
            
            if up_orf_test in orf5 and (down_orf_test in orf3 or down_orf_test1 in orf3):
                match5 = re.search(up_orf_test, orf5)
                aa5 = match5.end()
                                    
		if down_orf_test in orf3:
		    match3 = re.search(down_orf_test, orf3)
		    aa3 = match3.start() + 1
		else:
		    match3 = re.search(down_orf_test1, orf3)
		    aa3 = match3.start() + 1
                
                return (aa5, aa3, novel_aa)
                
    return False

def nuc_to_aa(nuc):
    if len(nuc) == 3:
	seq = Seq(nuc)
	return seq[0:].translate()
    return None
