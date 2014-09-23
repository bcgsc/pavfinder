import re
import sys
from intspan import intspan
from optparse import OptionParser

class Alignment:
    def __init__(self, query=None, qstart=None, qend=None, 
                 target=None, tstart=None, tend=None,
                 strand=None):
        self.query = query
        self.qstart, self.qend = qstart, qend
        self.target = target
        self.tstart, self.tend = tstart, tend
        self.strand = strand
        self.blocks = []
        self.query_blocks = []
        self.cigar = None
        # used in gapped_align.py for indel discovery
        self.cigarstring = ''
        self.identity = None
        self.score = None
	self.query_len = None
        
    @classmethod
    def from_alignedRead(cls, read, bam):
        target = bam.getrname(read.tid)
        align = cls(query=read.qname, target=target)
        align.strand = '-' if read.is_reverse else '+'
        align.query_len = get_query_len_from_cigar(read.cigar)
        align.blocks, align.query_blocks = cigar_to_blocks(read.cigar, read.pos + 1, align.strand)
        align.set_start_end_from_blocks()
        align.cigar = read.cigar
        align.cigarstring = read.cigarstring
	align.sam = read
        return align
    
    def set_start_end_from_blocks(self):
        if self.query_blocks and self.blocks:
            if self.query_blocks[0][0] < self.query_blocks[0][1]:
                self.qstart, self.qend = sorted([self.query_blocks[0][0], self.query_blocks[-1][1]], key=int)
                self.tstart, self.tend = sorted([self.blocks[0][0], self.blocks[-1][1]], key=int)
            else:
                self.qstart, self.qend = sorted([self.query_blocks[-1][1], self.query_blocks[0][0]], key=int)
                self.tstart, self.tend = sorted([self.blocks[-1][1], self.blocks[0][0]], key=int)
                
    def is_valid(self):
        try:
            assert self.qstart > 0
            assert self.qend > self.qstart
            assert self.tstart > 0
            assert self.tend > self.tstart
        except AssertionError:
            return False
        return True
    
    def as_bed(self):
	cols = [self.target, 
	        int(self.tstart) - 1,
	        int(self.tend),
	        self.query,
	        0,
	        '+',
	        int(self.tstart) - 1,
	        int(self.tend),
	        0,
	        len(self.blocks),
	        ','.join([str(int(b[1]) - int(b[0]) + 1) for b in self.blocks]),
	        ','.join([str(int(b[0]) - int(self.tstart)) for b in self.blocks])
	        ]
	        
	return '\t'.join([str(col) for col in cols])
        
def cigar_to_blocks(cigar, tstart, strand):
    query_len = get_query_len_from_cigar(cigar)
    qstart = 1 if strand == '+' else query_len
        
    tblocks = []
    qblocks = []          
    for i in range(len(cigar)):
        op, length = cigar[i]
        
        # 'D' (deletion)
        if op == 2 and i == 0:
            return None, None
        
        # 'S' or 'H' (clips)
        if op == 4 or op == 5:
            if i == 0:
                qstart = qstart + length if strand == '+' else qstart - length
                continue
        
        tblock = None
        qblock = None
        
        if not tblocks and op != 0:
            return None, None
        
        # match
        if op == 0:
            tend = tstart + length - 1
            qend = qstart + length - 1 if strand == '+' else qstart - length + 1
            tblock = [tstart, tend]
            qblock = [qstart, qend]
            
        # intron ('N'), skipped reference or deletion in reference ('D')
        elif op == 2 or op == 3:
            #intron
            if op == 3 and length < 3:
                continue
            
            qend = qend
            tend = tstart + length - 1 
            
        # insertion ('I') to reference
        elif op == 1:
            tend = tend
            qend = qstart + length - 1 if strand == '+' else qstart - length + 1
        
        if tblock:
            tblocks.append(tblock)
        if qblock:
            qblocks.append(qblock)
                    
        tstart = tend + 1
        qstart = qend + 1 if strand == '+' else qend - 1
                                     
    return tblocks, qblocks

def get_query_len_from_cigar(cigar):
    lens = []
    # if op is not 'D'(deletion), 'N'(skipped region), 'P'(padding)
    query_lens = [int(op[1]) for op in cigar if op[0] != 2 and op[0] != 3 and op[0] != 6]
    return sum(query_lens)

def get_match_len_from_cigar(cigar):
    query_lens = [int(op[1]) for op in cigar if op[0] <= 1]
    return sum(query_lens)

def reverse_complement(seq):
    """Reverse complements sequence string"""
    from string import maketrans
    complement = maketrans("agtcAGTC", "tcagTCAG")
    return seq[::-1].translate(complement)

def target_non_canonical(target):
    """Checks if target is non-canonical (haploptypes)"""
    if '_' in target:
	return True
    else:
	return False
    
def compare_chr(chr1, chr2):
    """For sorting chromosome names ignoring 'chr'"""
    if chr1[:3].lower() == 'chr':
	chr1 = chr1[3:]
    if chr2[:3].lower() == 'chr':
	chr2 = chr2[3:]
    
    if re.match('^\d+$', chr1) and not re.match('^\d+$', chr2):
	return -1
    
    elif not re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
	return 1
    
    else:
	if re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
	    chr1 = int(chr1)
	    chr2 = int(chr2)
	    
	if chr1 < chr2:
	    return -1
	elif chr1 > chr2:
	    return 1
	else:
	    return 0
	
def bam2track(in_file, out_file, track_name, desc=None):
    """Given a BAM file, output unique alignments in UCSC bed track
    
    Args:
        in_file: (Str) Path of BAM file
	out_file: (Str) Path of output file
	track_name: (Str) Name of track
	desc (optional): (Str) Description for track
	                 If not given, description will be same as name     
    
    Raises:
        IOError: when BAM file or output file is not valid
    """
    from itertools import groupby
    import pysam
    import sys
    
    try:
	bam = pysam.Samfile(in_file, 'rb')
    except IOError:
	sys.exit("can't open BAM:%s" % in_file)
	
    if desc is None:
	desc = track_name
	
    try:
	out = open(out_file, 'w')
    except IOError:
	sys.exit("can't open output file for writing:%s" % out_file)
	
    out.write('track name=%s description="%s"\n' % (track_name, desc))
    for query, group in groupby(bam.fetch(until_eof=True), lambda x: x.qname):
	alns = list(group)
	
	if len(alns) == 1 and not alns[0].is_unmapped:
	    align = Alignment.from_alignedRead(alns[0], bam)
	    try:
		out.write(align.as_bed() + '\n')
	    except:
		sys.stderr.write("can't convert alignment of contig %s to bed\n" % query)
    out.close()
    
if __name__ == '__main__':
    usage = "Usage: %prog"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-b", "--bam", dest="bam", help="alignment bam output")
    parser.add_option("-o", "--out_file", dest="out_file", help="output file")
    parser.add_option("-t", "--create_track", dest="create_track", help="create BED format UCSC track", action="store_true", default=False)
    parser.add_option("--name", dest="name", help="track name")
    parser.add_option("--desc", dest="desc", help="track description")
    
    (options, args) = parser.parse_args()
    if options.create_track:
        bam2track(options.bam, options.out_file, options.name, options.desc)