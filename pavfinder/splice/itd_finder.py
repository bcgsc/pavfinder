import re
import sys
import os
import subprocess
from pavfinder.shared.alignment import reverse_complement

class ITD_Finder: 
    @classmethod
    def detect_itd(cls, adj, align, contig_seq, outdir, min_len, max_apart, min_pid, debug=False):
	"""Determines if ins is an ITD
	
	Args:
	    adj: (Adjacency) original insertion adjacency
	    align: (Alignment) alignment
	    contig_seq: (str) contig sequence
	    outdir: (str) absolute path of output directory for temporarily storing blastn results if necessary
	    min_len: (int) minimum size of duplication
	    max_apart: (int) maximum distance apart of duplication
	    min_pid: (float) minimum percentage of identity of the 2 copies
	Returns:
	    If adj is determined to be an ITD, adj attributes will be updated accordingly
	"""
	novel_seq = adj.novel_seq if align.strand == '+' else reverse_complement(adj.novel_seq)
	dup = False
	# try regex first to see if perfect copies exist
	if len(novel_seq) >= min_len:
	    dup = cls.search_by_regex(novel_seq, contig_seq, max_apart)
	    if dup:
		cls.update_attrs(adj, align, dup)
	
	# try BLAST if regex fails
	if not dup:
	    dup = cls.search_by_align(adj.contigs[0], contig_seq, outdir, sorted(adj.contig_breaks[0]), min_len, max_apart, min_pid, debug=debug)
	    if dup:
		cls.update_attrs(adj, align, dup, contig_seq)
	
    @classmethod
    def update_attrs(cls, adj, align, dup, contig_seq=None):
	"""Updates attributes of original insertion adjacency to become ITD
	   1. rna_event -> 'ITD'
	   2. contig_breaks -> coordinates of the breakpoint (the last coordinate of the first and first coordinate of the second copy)
	   3. breaks -> genome coordinate of the last base of the first copy and the next position (if blast was done and <contig_seq>
	      given
	   4. novel_seq -> sequence of the second copy (if blast was done and <contig_seq> given
	   
	Args:
	    adj: (Adjacency) original insertion adjacency object
	    align: (Alignment) alignment obejct
	    dup: (tuple) a tuple of 2 tuples that contain the duplication coordinates in the contig
	    contig_seq(optional): (str) only given if break coordinates need to be modified (i.e. blastn has been done) and novel_seq
				  needs to be reset
	"""
	adj.rna_event = 'ITD'
	
	if contig_seq is not None:
	    # adjust genome break coordiantes based on new contig break coordinates
	    shift = dup[1][0] - adj.contig_breaks[0][0]
	    if align.strand == '+':
		new_genome_pos = adj.breaks[0] + shift
		new_genome_breaks = new_genome_pos - 1, new_genome_pos
	    else:
		new_genome_pos = adj.breaks[0] + shift * -1
		new_genome_breaks = new_genome_pos, new_genome_pos + 1
	    adj.breaks = new_genome_breaks
	    
	    # captures second copy of duplication as novel sequence
	    dup_size = max(dup[0][1] - dup[0][0] + 1, dup[1][1] - dup[1][0] + 1)
	    novel_seq = contig_seq[dup[1][0] - 1 : dup[1][0] - 1 + dup_size]
	    if align.strand == '-':
		novel_seq = reverse_complement(novel_seq)
	    adj.novel_seq = novel_seq

        
    @classmethod
    def search_by_regex(cls, novel_seq, contig_seq, max_apart):
	"""Finds if there is a duplication of the novel_seq within the contig sequence
	by simple regex matching i.e. only perfect matches captured
	
	Args:
	    novel_seq: (str) the original novel contig sequence defined in the insertion, assuming the 
			     length is already checked to be above the threshold
	    contig_seq: (str) the contig sequence
	    max_apart:  (int) the maximum number of bases allowed for the event to be called a duplication
	Returns:
	    a tuple of 2 tuples that contain the duplication coordinates in the contig
	"""
	matches = re.findall(novel_seq, contig_seq)
		
	if len(matches) > 1:
	    starts = []
	    p = re.compile(novel_seq)
	    for m in p.finditer(contig_seq):
		starts.append(m.start())
		
	    if len(starts) == 2 and\
	       starts[1] - (starts[0] + len(novel_seq)) <= max_apart:				
		return ((starts[0] + 1, starts[0] + len(novel_seq)),
		        (starts[1] + 1, starts[1] + len(novel_seq)))
	    
	return False
    
    @classmethod
    def parse_blast_tab(cls, blast_tab):
	"""Parses a tabular BLAST output file
	
	Args:
	    blast_tab: (str) full path to a tabular blastn output file (-outfmt 6 from blastn)
			     version of blastn: ncbi-blast-2.2.29+
	Returns:
	    A list of dictionary where the keys are the fields and the values are the results
	"""
	aligns = []
	fields = ('query', 'target', 'pid', 'alen', 'num_mismatch', 'num_gaps', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bit_score')
	for line in open(blast_tab, 'r'):
	    cols = line.rstrip('\n').split('\t')
	    # if somehow the number of values is different from the number of fields hard-coded above, no alignments will be captured
	    if len(cols) == len(fields):
		aligns.append(dict(zip(fields, cols)))
	    
	return aligns
	
    @classmethod
    def search_by_align(cls, contig, contig_seq, outdir, contig_breaks, min_len, max_apart, min_pid, debug=False):
	"""Checks if an insertion event is an ITD using blastn self-vs-self alignment
	Assumes:
	1. original contig coordiantes are the coordinates of the novel sequence
	2. original insertion will overlap the duplication ranges

	Conditions:
	1. percentage of identity between the duplications is at least <min_pid> (mismatch and gap allowed)
	2. the duplciations are at most <max_apart> apart
	3. the duplication is at least <min_len> long
	
	Returns:
	    a tuple of 2 tuples that contain the duplication coordinates in the contig
	"""						
	# run blastn of contig sequence against itself
	# generate FASTA file of contig sequence
	seq_file = '%s/tmp_seq.fa' % outdir
	out = open(seq_file, 'w')
	out.write('>%s\n%s\n' % (contig, contig_seq))
	out.close()
	
	# run blastn
	blast_output = '%s/tmp_blastn.tsv' % outdir
	cmd = 'blastn -query %s -subject %s -outfmt 6 -out %s' % (seq_file, seq_file, blast_output)

	try:
	    subprocess.call(cmd, shell=True)
	except:
	    # should check whether blastn is in the PATH right off the bet
	    sys.stderr.write('Failed to run BLAST:%s\n' % cmd)
	    sys.exit()
	    
	if os.path.exists(blast_output):
	    blast_aligns = cls.parse_blast_tab(blast_output)
	    # clean up
	    if not debug:
		for ff in (seq_file, blast_output):
		    os.remove(ff)
		
	    # identify 'significant' stretch of duplication from BLAST result
	    # and see if it overlaps contig break position
	    dups = []
	    # for picking the 'longest' dup if there are multiple candidates
	    dup_sizes = []
	    for aln in blast_aligns:
		# skip if a. it's the sequence match itself 
		#         b. match is too short 
		#         c. match is just a member of the 2 reciprocal matches
		if int(aln['alen']) == len(contig_seq) or\
		   int(aln['alen']) < min_len or\
		   int(aln['qstart']) >= int(aln['tstart']):
		    continue
		span1 = (int(aln['qstart']), int(aln['qend']))
		span2 = (int(aln['tstart']), int(aln['tend']))
		
		# 4 conditions for checking if it's an ITD
		#      1. longer than minimum length (checked above)
		#      2. percent of identity better than minumum
		#      3. the duplications less than or equal to max_apart
		#      4. either one the duplication ranges overlaps the original contigs breakpoint coordinate range
		#         (contig coordinates of the novel sequence)
		if float(aln['pid']) >= min_pid and\
		   abs(min(span2[1], span1[1]) - max(span2[0], span1[0])) <= max_apart and\
		   (min(contig_breaks[1], span1[1]) - max(contig_breaks[0], span1[0]) > 0 or
		    min(contig_breaks[1], span2[1]) - max(contig_breaks[0], span2[0]) > 0):
		    dup = [(int(aln['qstart']), int(aln['qend'])), 
		           (int(aln['tstart']), int(aln['tend']))]
		    dups.append(dup)
		    dup_sizes.append(int(aln['alen']))
		    
	    if dups:
		# if there are multiple candidates, pick the largest
		index_of_largest = sorted(range(len(dups)), key=lambda idx:dup_sizes[idx], reverse=True)[0]		
		dup = dups[index_of_largest]	
				
		return dup
		
	return False
