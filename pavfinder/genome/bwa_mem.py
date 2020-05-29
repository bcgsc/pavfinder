import re
from alignment import Alignment, reverse_complement, target_non_canonical
from subprocess import check_call, CalledProcessError
import sys
from intspan import intspan

def find_chimera(alns, bam, debug=False, check_haplotype=True):
    """Determine if given alignments are chimeric

    Args:
        alns: (List) List of Pysam AlignedRead objects
        bam: (AlignmentFile) Pysam handle to BAM file - for getting reference info
        debug: (Boolean) debug mode - will output debugging statements
        check_haplotype: (Boolean) whether to screen out alignments to references
                                   containing '_'
    """
    primary_alns = []
    secondary_alns = []
    for aln in alns:
        if re.search('[HS]', aln.cigarstring) and not aln.is_secondary:
            primary_alns.append(aln)
        else:
            secondary_alns.append(aln)
    
    if check_haplotype and len(primary_alns) > 1:
        replace_haplotype(primary_alns, secondary_alns, bam)
        
    if len(primary_alns) > 1:
        aligns = [Alignment.from_alignedRead(aln, bam) for aln in primary_alns]
        bad_aligns = [align for align in aligns if not align.is_valid()]
        if bad_aligns:
            if debug:
                for align in bad_aligns:
                    sys.stdout.write('bad alignment %s %s %s %s %s %s' % (align.query,
                                                                          align.qstart,
                                                                          align.qend,
                                                                          align.target,
                                                                          align.tstart,
                                                                          align.tend))
        else:
            valid_secondary_aligns = []
            if secondary_alns:
                secondary_aligns = [Alignment.from_alignedRead(aln, bam) for aln in secondary_alns]
                valid_secondary_aligns = [align for align in secondary_aligns if align.is_valid()]
                
            return aligns, valid_secondary_aligns
        
    return None, None
    
def replace_haplotype(primary_alns, secondary_alns, bam):
    """Replace primary alignments that maps to a haplotype contig with secondary alignments that don't
    
    Args:
        primary_alns: (list) Pysam AlignedRead objects of primary alignments
        alns: (list) Pysam AlignedRead objects of all alignments (primary and secondary)
        bam: Pysam bam handle for extracting the chromsome name
    """
    # indices of haplotype primary alignments to remove
    remove = set()
    for i in range(len(primary_alns)):
        chrom = bam.getrname(primary_alns[i].tid)
        # assume haplotypes have '_'
        if target_non_canonical(chrom):
            remove.add(i)
            
    # make changes to primary_alns only if it contains haplotype alignments
    if remove:
        # move haplotype alignments to secondary
        for i in sorted(remove, reverse=True):
            secondary_alns.append(primary_alns[i])
            del primary_alns[i]
            
        # move non-haplotype secondary alignments to primary
        for aln in secondary_alns:
            chrom = bam.getrname(aln.tid)
            # add secondary alignments that do not map to haplotypes for potential pairing
            if not target_non_canonical(chrom):
                primary_alns.append(aln)
                secondary_alns.remove(aln)
                
def effective_edit_distance(aln):
    """Returns edit distance without counting indels"""
    try:
        edit_distance = aln.opt('NM')
        indel_len = 0
        for op, length in aln.cigar:
            if op >= 1 and op <= 2:
                indel_len += length
                        
        return edit_distance - indel_len
    
    except:
        return None
            
def find_single_unique(alns, bam, debug=False):
    """Extracts single unique alignment for indel detection
    If there is only one alignment reported by BWA-mem even when '-a' is turned on
    
    Args:
        alns: (list) Pysam AlignedRead objects of the same contig
        bam: Pysam bam handle
    Returns:
        Alignment object or None
    """
    primary_alns = [aln for aln in alns if not aln.is_unmapped and not aln.is_secondary]
    if len(primary_alns) == 1:
        if primary_alns[0].mapq > 0:            
            matched_and_insertion_len = sum([a[1] for a in primary_alns[0].cigar if a[0] <= 1])
            if float(matched_and_insertion_len) / float(primary_alns[0].rlen) < 0.95:
                if debug:
                    sys.stdout.write('best alignment less than 0.95 mapped:%s %s\n' % (alns[0].qname, alns[0].cigarstring))
                return None
        
            else:
                edit_distance = effective_edit_distance(alns[0])
                if edit_distance is not None and float(edit_distance)/float(primary_alns[0].inferred_length) > 0.1:
                    if debug:
                        sys.stdout.write('filter out single uniq alignment %s: edit distance %s - > 0.1 of contig len %d (%.01f)\n' % (alns[0].qname,
                                                                                                                                       edit_distance,
                                                                                                                                       primary_alns[0].inferred_length,
                                                                                                                                       float(edit_distance)/float(primary_alns[0].inferred_length)
                                                                                                                                       ))
                    return None
                        
        else:
            if debug:
                sys.stdout.write('filter out single uniq alignment %s: mapq = 0\n' % primary_alns[0].qname)
            return None
            
        #ambiguous_NM = 5
        #for aln in alns:
            #if aln.is_secondary and \
               #not re.search('[HS]', aln.cigarstring) and\
               #re.match('\d+M', aln.cigarstring) and re.search('\d+M$', aln.cigarstring) and\
               #int(aln.opt('NM')) - int(primary_alns[0].opt('NM')) <= ambiguous_NM:
                #if debug:
                    #sys.stdout.write('secondary alignments too similar %s\n' % primary_alns[0].qname)
                #return None
        
        return Alignment.from_alignedRead(primary_alns[0], bam) 
    else:
        return None

def run(fasta, output_bam, genome, index_dir=None, num_threads=4):
    """Runs BWA-mem on command line
    
    Args:
        fasta: (str) Path of input Fasta file
        output_bam: (str) Path of output BAM
        genome: (str) Prefix of genome indices in index_dir
        index_dir: (str) Path of directory containing indices of genome
        num_threads: (int) Number of threads in running alignment
    Returns:
        return code of system call
    """
    cmd = 'bwa mem -a -t %d %s/%s %s | samtools view -bhS - -o %s' % (num_threads, index_dir, genome, fasta, output_bam)

    print(cmd)
    try:
        returncode = check_call(cmd, shell=True)
    except CalledProcessError, e:
        sys.stderr.write('Failed to align:%s\n' % e.cmd)
        returncode = e.returncode
        
    return returncode

def find_microhomology(aligns, contig_seq):
    """Finds micromology given 2 alignments and contig sequence
    The homology sequence is based on the contig sequence
    
    Homology is found based on the fact that BWA-mem will report overlapping
    contig coordinates in chimeric alignments.
    
    Args:
        aligns: (list) 2 Alignment objects of chimera
        contig_seq: (str) Contig sequence
    Returns:
        Tuple of homology sequence(str) and homology (contig)coordinates((int, int)) 
    """
    homol_seq = None
    homol_coords = None
    
    contig_span1 = intspan('%s-%s' % (aligns[0].qstart, aligns[0].qend))
    contig_span2 = intspan('%s-%s' % (aligns[1].qstart, aligns[1].qend))
    overlap = contig_span1.intersection(contig_span2)
    if len(overlap) > 0:
        homol_coords = overlap.ranges()[0]
        homol_seq = contig_seq[homol_coords[0] - 1 : homol_coords[1]]
    
    return homol_seq, homol_coords

def find_untemplated_sequence(aligns, contig_seq):
    """Finds untemplated sequence in chimeric breakpoint
    This corresponds to any sequence at the breakpoint that is not covered by the
    2 alignments in the chimera.  
    The sequence will be given is the same strand as the first alignment 
    (the first and second alignments should have the same strand)
    
    Args:
        aligns: (list) 2 Alignment objects of chimera
        contig_seq: (str) Contig sequence
    Returns:
        Untemplated sequence or None
    """
    untemplated_seq = '-'
    
    contig_span1 = intspan('%s-%s' % (aligns[0].qstart, aligns[0].qend))
    contig_span2 = intspan('%s-%s' % (aligns[1].qstart, aligns[1].qend))
    sorted_contig_coords = sorted([aligns[0].qstart, aligns[0].qend, aligns[1].qstart, aligns[1].qend])
    whole_span = intspan('%s-%s' % (min(sorted_contig_coords), max(sorted_contig_coords)))
    unmapped = whole_span - contig_span1 - contig_span2
    
    if len(unmapped) > 0:
        unmapped_coords = unmapped.ranges()
        untemplated_seq = contig_seq[unmapped_coords[0][0] - 1 : unmapped_coords[0][1]]
        # sequence given in relation to strand of first alignment
        if aligns[0].strand == '-':
            untemplated_seq = reverse_complement(untemplated_seq)
    
    return untemplated_seq
