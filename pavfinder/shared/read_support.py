import pysam
import sys
import re
from itertools import groupby, chain
from optparse import OptionParser
from sets import Set
import multiprocessing as mp
from intspan import intspan
from collections import OrderedDict, defaultdict
from operator import itemgetter

def find_flanking(reads, span, contig_len, overlap_buffer=1, allow_clipped=False, min_ratio_mapped=None, debug=False):
    uniq_frags = defaultdict(list)

    tlens = []
    proper_pairs = Set()
    for read in reads:
	if read.is_proper_pair and is_fully_mapped(read, contig_len, allow_clipped=allow_clipped, min_ratio_mapped=min_ratio_mapped):
	    proper_pairs.add(read.qname + str(read.pos))
	    #if read.tlen > 0:
		#tlens.append(read.tlen)

    for read in [r for r in reads if not r.is_unmapped and r.tlen > 0]:
	if read.qname + str(read.pos) in proper_pairs:
	     # internal fragment
	    frag = (read.pos + read.alen, read.pnext)

	    if frag is not None and frag[0] <= span[0] - overlap_buffer and frag[1] >= span[1] + overlap_buffer:
		uniq_frags[frag].append(read.qname)

    num_frags = len(uniq_frags.keys())

    if debug:
	for f in uniq_frags.keys():
	    sys.stdout.write("Accepted flanking: %s %s %s\n" % (span, f, ','.join(uniq_frags[f])))

    return num_frags, tlens

def find_spanning(reads, span, contig_seq, overlap_buffer=1, debug=False, perfect=False, get_seq=False, allow_clipped=False, min_ratio_mapped=None):
    def check_read(start, end):
	break_seq = contig_seq[start - 1 - overlap_buffer: end + overlap_buffer]
	return read.pos < start - overlap_buffer and\
	       read.pos + read.alen >= end + overlap_buffer and\
	       is_fully_mapped(read, contig_len, perfect=perfect, allow_clipped=allow_clipped, min_ratio_mapped=min_ratio_mapped) and\
	       is_break_region_perfect(read, break_seq, (start, end), overlap_buffer)

    contig_len = len(contig_seq)
    pos = Set()
    names = Set()
    for read in reads:
	if read.alen:
	    if check_read(span[0], span[1]):
		strand = '+' if not read.is_reverse else '-'

		start_pos = read.pos
		if read.cigar[0][0] >= 4 and read.cigar[0][0] <= 5:
		    start_pos = -1 * read.cigar[0][1]
		key = str(start_pos) + strand
		
		if not key in pos and not read.qname in names:
		    pos.add(key)
		    names.add(read.qname)

		    if debug:
			sys.stdout.write("Accepted spanning read(perfect:%s): %s %s %s %s %s\n" % (perfect, 
		                                                                                   read.qname,
		                                                                                   span,
		                                                                                   (read.pos + 1, read.pos + read.alen),
		                                                                                   read.seq,
		                                                                                   strand))
    num_frags = len(pos)

    return num_frags

def is_break_region_perfect(read, break_seq, breaks, overlap_buffer):
    start_idx = breaks[0] - read.pos - 1
    if read.cigar[0][0] >= 4 and read.cigar[0][0] <= 5:
	start_idx += read.cigar[0][1]
    read_break_seq = read.seq[start_idx - overlap_buffer : start_idx - overlap_buffer + len(break_seq)]
    if read_break_seq.lower() == break_seq.lower():
	return True
    else:
	return False
    
def check_tiling(reads, breaks, contig_len, debug=False):
    """Checks if there are reads tiling across breakpoints with no gaps
    
    This will be used for checking integrity of breakpoint where there is a novel
    sequence of considerable size and there is not enough flanking sequences for read pairs
    to suggest validity of fragment
    
    Args:
        reads:  (list) Pysam AlignedRead objects
        breaks: (tuple) sorted coordinates of breakpoint positions in contigs
    Returns:
        Boolean if there are reads spanning across breakpoints with no gaps
    """
    span = None
    for read in reads:
	# skip reads that is unmapped, not properly paired, or the second mate, or not fully mapped
        if not read.alen or not is_fully_mapped(read, contig_len):
            continue
	
	# skip reads that don't overlap the breakpoints
	if read.pos + read.alen < breaks[0] or read.pos > breaks[1]:
	    continue
	
	try:
	    span = span.union(intspan('%d-%d' % (read.pos + 1, read.pos + read.alen)))
	except:
	    span = intspan('%d-%d' % (read.pos + 1, read.pos + read.alen))
	    	
    if span is not None:
	break_span = intspan('%d-%d' % (breaks[0], breaks[1]))
	# make sure there is no gap in tiling reads and spans the entire breakpoint
	if len(span.ranges()) == 1 and len(span & break_span) == len(break_span):
	    return True
    
    return False

def is_fully_mapped(read, contig_len, perfect=False, allow_clipped=False, min_ratio_mapped=None):
    if not re.search('[HS]', read.cigarstring) or is_fully_mapped_to_edge(read, contig_len):
        if perfect:
            return True if read.opt('NM') == 0 else False
        return True
    
    elif allow_clipped and\
         min_ratio_mapped is not None and\
         len(read.cigar) == 2 and\
         re.search('S', read.cigarstring) and\
         float(read.alen) / read.inferred_length >= min_ratio_mapped:
	if perfect:
            return True if read.opt('NM') == 0 else False
	return True

    return False

def is_fully_mapped_to_edge(read, contig_len):
    if read.cigar:
	if len(read.cigar) == 2:
	    # clipped at beginning
	    if (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) and read.cigar[1][0] == 0 and read.pos == 0 and read.inferred_length == read.rlen:
		return True
	    # clipped at the end
	    elif read.cigar[0][0] == 0 and (read.cigar[1][0] or 4 and read.cigar[1][0] == 5) and read.pos + read.alen == contig_len and read.inferred_length == read.rlen:
		return True
	elif len(read.cigar) == 3 and\
	     read.pos == 0 and read.pos + read.alen == contig_len and\
	     (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) and\
	     read.cigar[1][0] == 0 and\
	     (read.cigar[2][0] == 4 or read.cigar[2][0] == 5) and\
	     read.inferred_length == read.rlen:
	    return True
    
    return False

def worker(args):
    bam_file, coords, overlap_buffer, contig_fasta_file, perfect, get_seq, allow_clipped, min_ratio_mapped, debug = args
    contig_fasta = pysam.Fastafile(contig_fasta_file)

    coords_batch = defaultdict(list)
    for contig, start, end, align_type in coords:
	coords_batch[contig].append(((start, end), align_type))
    
    results, tlens = fetch_support(coords_batch,
                                   bam_file,
                                   contig_fasta,
                                   overlap_buffer=overlap_buffer,
                                   perfect=perfect,
                                   get_seq=get_seq,
                                   allow_clipped=allow_clipped,
                                   min_ratio_mapped=min_ratio_mapped,
                                   debug=debug)

    supports = []
    for contig, coords in results.iteritems():
	for coord in coords:
	    start, end = map(int, coord.split('-'))
	    supports.append([contig, start, end,
	                     results[contig][coord][0],
	                     results[contig][coord][1],
	                     results[contig][coord][2],
	                     results[contig][coord][3]])
    supports.append(tlens)
    return supports

def create_batches(bam_file, coords, size, overlap_buffer, contig_fasta_file, perfect, get_seq, allow_clipped, min_ratio_mapped, debug=False):
    if size == 0:
        yield bam_file, coords, overlap_buffer, contig_fasta_file, perfect, get_seq, debug
    else:
        for i in xrange(0, len(coords), size):
            if len(coords) - (i + size) < size:
                yield bam_file, coords[i:], overlap_buffer, contig_fasta_file, perfect, get_seq, allow_clipped, min_ratio_mapped, debug
                break
            else:
                yield bam_file, coords[i:i + size], overlap_buffer, contig_fasta_file, perfect, get_seq, allow_clipped, min_ratio_mapped, debug

def fetch_support(coords, bam_file, contig_fasta, overlap_buffer=0, perfect=False, get_seq=False, allow_clipped=False, min_ratio_mapped=None, debug=False):
    """Fetches read support when number given coords is relatively small
    It will use Pysam's fetch() instead of going through all read alignments
    Args:
        coords: (dictionary) coords[contig] = [spans]
                spans = (list) of (start, end) where 'start' and 'end' are not sorted
        bam_file: (string) absolute path of reads-to-contigs bam file
        debug: (boolean) prints debug statements
    Returns:
        dictionary of number of read support
        results[contig][coords] = (spanning_pos, spanning_neg, flanking)
        contig: (str)contig name
        coords: (str) original given coordinates (won't be sorted) "coord1-coord2"
        spanning_pos: (int) number unique postively spanning reads
        spanning_neg: (int) number unique postively spanning reads
        flanking: (int) number of unique flanking pairs
    """
    bam = pysam.Samfile(bam_file, 'rb')
    results = {}
    tlens_all = []
    count = 0
    window_size = 2000
    for contig, spans_align_types in coords.iteritems():
	count += 1
        results[contig] = {}
	contig_seq = contig_fasta.fetch(contig)
	contig_len = len(contig_seq)

	spans = []
	for span, align_type in spans_align_types:
	    spans.append(span)
	    coord = '-'.join(map(str, span))
	    results[contig][coord] = [0, 0, None, []]

	    reads = []
	    for read in bam.fetch(contig, max(0, span[0] - window_size), min(contig_len, span[1] + window_size)):
		reads.append(read)

	    spanning = find_spanning(reads, span, contig_seq, debug=debug,
	                             overlap_buffer=overlap_buffer,
	                             perfect=perfect,
	                             get_seq=get_seq,
	                             allow_clipped=allow_clipped,
	                             min_ratio_mapped=min_ratio_mapped)

	    flanking, tlens = find_flanking(reads,
	                                    span,
	                                    contig_len,
	                                    overlap_buffer=overlap_buffer,
	                                    allow_clipped=allow_clipped,
	                                    min_ratio_mapped=min_ratio_mapped,
	                                    debug=debug)

	    results[contig][coord][0] = spanning
	    results[contig][coord][1] = flanking
	#print count, len(coords.keys()), contig_len

    return results, tlens_all

def expand_contig_breaks(chrom, breaks, contig, contig_breaks, event, ref_fasta, contig_fasta, debug=False):
    def extract_repeat(seq):
	repeat = {'start':None, 'end':None}
	if len(seq) == 1:
	    repeat['start'] = seq
	    repeat['end'] = seq
	else:
	    re_start = re.compile(r"^(.+?)\1+")
	    re_end = re.compile(r"(.+?)\1+$")

	    repeats = re_start.findall(seq)
	    if repeats:
		repeat['start'] = repeats[0]
	    repeats = re_end.findall(seq)
	    if repeats:
		repeat['end'] = repeats[0]
	
	return repeat

    contig_breaks_sorted = sorted(contig_breaks)
    pos_strand = True if contig_breaks[0] < contig_breaks[1] else False
    contig_seq = contig_fasta.fetch(contig)
    contig_breaks_expanded = [contig_breaks[0], contig_breaks[1]]
    
    seq = None
    if event == 'del':
	seq = ref_fasta.fetch(chrom, breaks[0], breaks[1] - 1)
	if not pos_strand:
	    seq = reverse_complement(seq)
	    
    elif event == 'ins':
	seq = contig_seq[contig_breaks_sorted[0] : contig_breaks_sorted[1] - 1]
	    
    if seq is None:
	return None
    
    repeat = extract_repeat(seq)
    
    if repeat['end'] is not None:
	seq = repeat['end']
	size = len(seq)
	# downstream
	start = contig_breaks_sorted[1]
	expand = 0
	while start + size <= len(contig_seq):
	    next_seq = contig_seq[start : start + size]
	    if next_seq.upper() != seq.upper():
		break
	    else:
		expand += size
		start += size
	if pos_strand:
	    contig_breaks_expanded[1] += expand
	else:
	    contig_breaks_expanded[0] += expand
	
    # upstream
    if repeat['start'] is not None:
	seq = repeat['start']
	size = len(seq)
	start = contig_breaks_sorted[0] - 1
	expand = 0
	while start - size >= 0:
	    next_seq = contig_seq[start - size : start]
	    if next_seq.upper() != seq.upper():
		break
	    else:
		expand -= size
		start -= size
	if pos_strand:
	    contig_breaks_expanded[0] += expand
	else:
	    contig_breaks_expanded[1] += expand
	
    if debug and tuple(contig_breaks_expanded) != contig_breaks:
	sys.stdout.write('contig breaks expanded:%s %s -> %s\n' % (contig, contig_breaks, contig_breaks_expanded))
	
    return tuple(contig_breaks_expanded)

def sum_support(values, use_minimum=False):
    if not use_minimum:
	return sum(values)
    else:
	return min(values)

def filter_support(spanning, flanking, min_support, use_spanning=True, use_flanking=False):
    support = 0
    if use_spanning and use_flanking:
	support = spanning + flanking
    elif use_spanning:
	support = spanning
    elif use_flanking:
	support = flanking

    return support >= min_support

def main(args, options):
    coords_file = args[0]
        
    coords = []
    for line in open(coords_file, 'r'):
        cols = line.rstrip('\n').split()
        span = tuple(sorted((int(cols[1]), int(cols[2]))))
        #coords[cols[0]].append((span, cols[3]))
	coords.append((cols[0], span[0], span[1], cols[3]))
    coords_sorted = sorted(coords, key=itemgetter(0,1,2))
            
    batches = list(create_batches(args[1],
                                  coords_sorted,
                                  len(coords_sorted)/options.num_procs,
                                  options.min_overlap,
                                  args[2],
                                  False,
                                  False,
                                  options.allow_clipped_support,
                                  options.support_min_mapped,
                                  ))
    print 'uuu', len(batches)
    pool = mp.Pool(processes=options.num_procs)
    batch_results = pool.map(worker, batches)
    pool.close()
    pool.join()
    
    results = {}
    tlens_all = []
    support_reads = []
    for batch_result in batch_results:
        for support in batch_result:
	    # insert sizes
	    if len(support) > 7:
		tlens_all.extend(support)
		continue

	    if len(support) == 0:
		continue

            contig, start, stop, spanning, flanking, tiling, support_reads = support
            coords = '%s-%s' % (start, stop)
            try:
                results[contig][coords] = (spanning, flanking, tiling, support_reads)
            except:
                results[contig] = {}
                results[contig][coords] = (spanning, flanking, tiling, support_reads)

    with open(args[3], 'w') as out:
	for contig in results.keys():
	    for coords in results[contig]:
		out.write('%s\n' % '\t'.join(map(str, [contig, coords, results[contig][coords][0], results[contig][coords][1]])))

if __name__ == '__main__':
    usage = "Usage: %prog coords_file bamfile contigs_fasta out_file"
    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--num_procs", dest="num_procs", help="Number of processes. Default: 5", default=5, type=int)
    parser.add_option("--min_overlap", dest="min_overlap", help="minimum breakpoint overlap for identifying read support. Default:4", type='int', default=4)
    parser.add_option("--allow_clipped_support", dest="allow_clipped_support", help="allow using clipped reads in gathering read support", action="store_true", default=False)
    parser.add_option("--support_min_mapped", dest="support_min_mapped", help="when clipped reads are allowed as read support, minimum ratio of read length mapped Default:0.8", type='float', default=0.8)
    (options, args) = parser.parse_args()
    if len(args) == 4:
        main(args, options)