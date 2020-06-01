import re
import pysam
from .alignment import Alignment, reverse_complement
from .adjacency import Adjacency
from . import bwa_mem


def find_single_unique(alns, bam, debug=False):
    return bwa_mem.find_single_unique(alns, bam, debug=debug)


def find_adjs(align, contig_seq, is_transcriptome, ins_as_ins=False, query_fasta=None, target_fasta=None):
    adjs = []

    if not align.is_valid() or not align.cigarstring or len(align.blocks) != len(align.query_blocks) or\
       not re.search(r'[ID]\d+', align.cigarstring):
        return adjs

    # identify all (target) gaps in cigar: N=intron, D=deletion
    count = 0
    target_gaps = {}
    for gap in re.findall(r'[ND]\d+', align.cigarstring):
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
            if gap_count in target_gaps:
                breaks = (align.blocks[i][1], align.blocks[i + 1][0])
                contig_breaks = (align.query_blocks[i][1], align.query_blocks[i + 1][0])

                del_seq = target_fasta.fetch(align.target, breaks[0], breaks[1] - 1)
                if align.strand == '-':
                    del_seq = reverse_complement(del_seq)

                duplicated, repeat_seq, repeat_num = is_duplicated(del_seq, sorted(contig_breaks, key=int), contig_seq)
                if duplicated[0] > 0 or duplicated[1] > 0:
                    contig_breaks = (contig_breaks[0] - len(repeat_seq) * duplicated[0],
                                     contig_breaks[1] + len(repeat_seq) * duplicated[1])

                if target_gaps[gap_count][0] == 'N' and is_transcriptome:
                    continue

                if query_gap == 0:
                    rearrangement = 'del'
                else:
                    rearrangement = 'indel'

                probe_seq, break_pos = Adjacency.extract_probe(contig_seq, contig_breaks)

                adj = Adjacency((align.target, align.target),
                                breaks,
                                rearrangement,
                                contig=align.query,
                                contig_sizes=len(contig_seq),
                                contig_breaks=contig_breaks,
                                probes=probe_seq,
                                orients=('L', 'R'),
                                aligns=[align],
                                align_types='gapped')
                adjs.append(adj)

                if duplicated[0] > 0 or duplicated[1] > 0:
                    if align.strand == '+':
                        adj.repeat_seq = repeat_seq.upper()
                    else:
                        adj.repeat_seq = reverse_complement(repeat_seq.upper())
                    adj.repeat_num = repeat_num
                    adj.repeat_num_change = '{}>{}'.format(duplicated[0] + duplicated[1] + repeat_num,
                                                           duplicated[0] + duplicated[1])

            gap_count += 1

        # insertion
        elif query_gap > 0:
            event = 'ins'
            breaks = (align.blocks[i][1], align.blocks[i][1])
            contig_breaks = (align.query_blocks[i][1], align.query_blocks[i + 1][0])

            if align.strand == '+':
                novel_seq = contig_seq[align.query_blocks[i][1]: align.query_blocks[i][1] + query_gap]
                novel_seq_ref = novel_seq
            else:
                novel_seq = contig_seq[align.query_blocks[i + 1][0]: align.query_blocks[i + 1][0] + query_gap]
                novel_seq_ref = reverse_complement(novel_seq)

            length_novel_seq = len(novel_seq)

            duplicated, repeat_seq, repeat_num = is_duplicated(novel_seq, sorted(contig_breaks, key=int), contig_seq)
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

                contig_breaks = (contig_breaks[0] - len(repeat_seq) * duplicated[0],
                                 contig_breaks[1] + len(repeat_seq) * duplicated[1])

            probe_seq, break_pos = Adjacency.extract_probe(contig_seq, contig_breaks)

            adj = Adjacency((align.target, align.target),
                            breaks,
                            event,
                            contig=align.query,
                            contig_breaks=contig_breaks,
                            contig_sizes=len(contig_seq),
                            probes=probe_seq,
                            novel_seq=novel_seq_ref,orients=('L', 'R'),
                            aligns=[align],
                            align_types='gapped')

            if event == 'dup':
                if align.strand == '+':
                    adj.repeat_seq = repeat_seq.upper()
                else:
                    adj.repeat_seq = reverse_complement(repeat_seq.upper())
                adj.repeat_num = repeat_num
                adj.repeat_num_change = '{}>{}'.format(duplicated[0] + duplicated[1], 
                                                       duplicated[0] + duplicated[1] + repeat_num)

            adjs.append(adj)

    return adjs


def is_homopolymer(seq):
    return len(set(map(''.join, zip(*[iter(seq)] * 1)))) == 1


def find_repeat(seq):
    def chop_seq(seq, size):
        """chop seq into equal size sub-strings"""
        return map(''.join, zip(*[iter(seq)] * size))

    size_range = [2, 4]
    for size in range(size_range[0], size_range[1] + 1):
        if len(seq) % size == 0:
            uniq_subseqs = set(chop_seq(seq, size))
            if len(uniq_subseqs) == 1:
                return list(uniq_subseqs)[0]

    return None


def is_duplicated(novel_seq, contig_breaks, contig_seq, min_len=3):
    repeat_seq, repeat_num = None, None
    if is_homopolymer(novel_seq):
        repeat_seq = novel_seq[0]
        repeat_num = len(novel_seq)
    else:
        repeat_seq = find_repeat(novel_seq)
        if repeat_seq is not None:
            repeat_num = int(len(novel_seq) / len(repeat_seq))

    duplicated = [0, 0]
    if repeat_seq is not None:
        repeat_len = len(repeat_seq)

        # upstream
        num_repeats = 0
        start = contig_breaks[0] - repeat_len
        while start >= 0:
            before_seq = contig_seq[start: start + repeat_len]
            if repeat_seq.upper() == before_seq.upper():
                num_repeats += 1
                start -= repeat_len
            else:
                break
        duplicated[0] = num_repeats

        # downstream
        num_repeats = 0
        start = contig_breaks[1] - 1
        while start + repeat_len < len(contig_seq):
            after_seq = contig_seq[start: start + repeat_len]
            if repeat_seq.upper() == after_seq.upper():
                num_repeats += 1
                start += repeat_len
            else:
                break
        duplicated[1] = num_repeats

    return duplicated, repeat_seq, repeat_num


def screen_probe_alns(adj_aligns, probe_alns, align_type, min_pc_mapped=1.0):
    for aln in probe_alns:
        if aln.is_unmapped:
            continue

        matched_len = 0
        query_len = sum([a[1] for a in aln.cigar if a[0] in (0, 1, 4, 5)])
        matched_len = sum([a[1] for a in aln.cigar if a[0] == 0])
        # allows gap inside probe's alignment for checking end-to-end mapping
        if re.match(r'^\d+M.+\d+M$', aln.cigarstring):
            matched_len = query_len

        # if align_type == 'split' and (re.match('\d+M$', aln.cigarstring) or
        # re.match('\d+M.+\d+M$', aln.cigarstring)):
        if align_type == 'split' and float(matched_len) / float(query_len) >= min_pc_mapped:
            return False

        # if it's a single alignment and the probe can map perfectly (no clips,
        # no insertion/deletion) to a location -> out
        if align_type == 'gapped' and aln.rlen == aln.alen and not re.search('[DIN]', aln.cigarstring):
            return False

    return True
