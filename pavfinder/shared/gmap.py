from alignment import Alignment
from subprocess import check_call, CalledProcessError
import sys

def find_chimera(alns, bam, debug=False):
    num_chimeric = 0
    
    for aln in alns:
        chimeric = [field for field in aln.tags if field[0] == 'XT']
        if chimeric and len(chimeric) == 1:
            num_chimeric += 1
          
    if num_chimeric == len(alns):
        aligns = [Alignment.from_alignedRead(aln, bam) for aln in alns]
        
        # check if all alignments are valid
        for align in aligns:
            if not align.is_valid():
                print 'problem', align.query, align.qstart, align.qend, align.tstart, align.tend
                return []

        return aligns
    else:
        return []
    
def find_single_unique(alns, bam, debug=False):
    if len(alns) == 1 and not alns[0].is_unmapped:
        return Alignment.from_alignedRead(alns[0], bam) 
    
    return None
        
def is_exon_exon_fusion(alns, prob_cutoff=0.9):
    num_high_prob = 0
    for aln in alns:
        chimeric = [field for field in aln.tags if field[0] == tag]
        splice_sites, prob_donor, prob_acceptor, lens = chimeric[0][1].split(',')
        if float(prob_donor) >= prob_cutoff and float(prob_acceptor) >= prob_cutoff:
            num_high_prob += 1
               
    return num_high_prob == len(alns)

def run(fasta, output_bam, genome, index_dir=None, num_threads=4, multi=False):
    if index_dir is None:
        if not multi:
            cmd = 'gmap -d %s %s -t %d -A -f samse -O -x 10 -n 0 | samtools view -bhS - -o %s' %\
                (genome, fasta, num_threads, output_bam)
        else:
            cmd = 'gmap -d %s %s -t %d -A -f samse -O -x 10 | samtools view -bhS - -o %s' %\
                (genome, fasta, num_threads, output_bam)
    else:
        if not multi:
            cmd = 'gmap -d %s -D %s %s -t %d -A -f samse -O -x 10 -n 0 | samtools view -bhS - -o %s' %\
                (genome, index_dir, fasta, num_threads, output_bam)
        else:
            cmd = 'gmap -d %s -D %s %s -t %d -A -f samse -O -x 10 | samtools view -bhS - -o %s' %\
                (genome, index_dir, fasta, num_threads, output_bam)

    print cmd
    try:
        returncode = check_call(cmd, shell=True)
    except CalledProcessError, e:
        sys.stderr.write('Failed to align:%s\n' % e.cmd)
        returncode = e.returncode
        
    return returncode

def find_microhomology(aln, contig_seq):
    homol_seq = None
    homol_coords = None
    
    try:
        end_base1, end_base2 = aln.opt('XT').split(',')[-1].split('..')
        if end_base1 != end_base2 and end_base1.isdigit() and end_base2.isdigit():
            homol_coords = (int(end_base1) + 1, int(end_base2))
            homol_seq = contig_seq[int(end_base1):int(end_base2)]
    except:
        print 'XT tag empty for contig %s' % aln.qname
        
    return homol_seq, homol_coords

def is_chimera(aln):
    return aln.has_tag('XT')