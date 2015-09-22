import re
import sys
import os
import subprocess
import pysam
from itertools import groupby
from pavfinder.shared.alignment import reverse_complement, Alignment
from pavfinder.splice.event import Event
import pavfinder.SV.gapped_align as gapped_align
import pavfinder.SV.split_align as split_align
from intspan import intspan

conditions = None
contigs_to_transcripts = None

def detect_itd(adj, align, contig_seq, outdir, min_len, max_apart, min_pid, debug=False):
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
	dup = search_by_regex(novel_seq, contig_seq, max_apart)
	if dup:
	    update_attrs(adj, align, dup)
    
    # try BLAST if regex fails
    if not dup:
	dup = search_by_align(adj.contigs[0], contig_seq, outdir, sorted(adj.contig_breaks[0]), min_len, max_apart, min_pid, debug=debug)
	if dup:
	    update_attrs(adj, align, dup, contig_seq)
	    
def call_from_chimera(adj, junc_matches1, junc_matches2, transcripts, contig_seq, ref_fasta, outdir=None, debug=False):
    """Call ITD from split alignments

    Will update adj with ITD attributes if it is affirmed

    Args:
        adj: (Adjacency) Adjacent object of chimera
	junc_matches1: (dict) key=transcript_id, value=list of tuple(exon_index, match_operator) e.g. {'NM_004119': [(10, '=<')]}
	junc_matches2: (dict) key=transcript_id, value=list of tuple(exon_index, match_operator) e.g. {'NM_004119': [(10, '=<')]}
	transcripts: (dict) key=transcript_name value=Transcript object
    """
    txt_id1 = junc_matches1.keys()[0]
    txt_id2 = junc_matches2.keys()[0]
    if junc_matches1[txt_id1] is not None and junc_matches2[txt_id2] is not None:
	exon_idx1 = junc_matches1[txt_id1][0][0]
	exon_idx2 = junc_matches2[txt_id2][0][0]
	# breaks in same exon
	if txt_id1 == txt_id2 and exon_idx1 == exon_idx2:
	    txt = transcripts[txt_id1]
	    adj = confirm_with_transcript_seq(adj.contigs[0], contig_seq, txt, ref_fasta, adj=adj, outdir=outdir, debug=debug)
	
def update_attrs(adj, align, dup, contig_seq=None):
    """Updates adj with ITD-related attributes
    
    Used only by detect_itd()
    Args:
        adj: (Adjacency) original insertion adjacency object
	align: (Alignment) alignment obejct
	dup: (tuple) a tuple of 2 tuples that contain the duplication coordinates in the contig
	contig_seq(optional): (str) only given if break coordinates need to be modified (i.e. blastn has been done) and novel_seq
                              needs to be reset
    """
    adj.rna_event = 'ITD'
    novel_seq = adj.novel_seq

    terminal_bases_changed = False
    if len(novel_seq) <= len(adj.novel_seq) and len(adj.novel_seq) - len(novel_seq) <= 2:
	start_pos = adj.novel_seq.find(novel_seq)
	end_pos = start_pos + len(novel_seq) - 1
	if start_pos == 0 or end_pos == len(adj.novel_seq) - 1:
	    terminal_bases_changed = True

    if not terminal_bases_changed:
	novel_seq = contig_seq[dup[1][0] - 1 : dup[1][0] - 1 + dup_size]

    dup_size = len(adj.novel_seq)

    adj.breaks = [adj.breaks[0] + dup_size, adj.breaks[0] + 1]
    adj.contig_breaks = [(adj.contig_breaks[0][1] - 1, adj.contig_breaks[0][1])]
    adj.size = dup_size
    adj.novel_seq = None
	
def search_by_regex(novel_seq, contig_seq, max_apart):
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

def parse_blast_tab(blast_tab):
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

def search_by_align(contig, contig_seq, outdir, contig_breaks, min_len, max_apart, min_pid, debug=False):
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
	blast_aligns = parse_blast_tab(blast_output)
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

def confirm_with_transcript_seq(contig, contig_seq, transcript, ref_fasta, outdir='.', adj=None, debug=False):
    txt_seq = transcript.get_sequence(ref_fasta)

    if contigs_to_transcripts is not None:
	return capture_from_transcript_align(contigs_to_transcripts, contig, contig_seq, transcript, txt_seq, outdir, adj, debug=debug)
    else:
	return None

def capture_from_transcript_align(bam, contig, contig_seq, transcript, txt_seq, outdir, adj_genome=None, debug=True):
    """Captures ITD from contig-transcript BWA-mem alignment

    Args:
        bam: (AlignmentFile) contig-transcript BWA-mem alignment Pysam handle
	contig: (str) contig ID
	contig_seq: (str) contig sequence
	transcript: (Transcript) Transcript object
	txt_seq: (str) transcript sequence
	adj_genome: (Adjacency) Adjacency from genome alignment
    """
    itds = []
    # grep alignments mapped to transcript and sort by query(contig) name
    txt_alns = []
    for aln in bam.fetch(transcript.id):
	txt_alns.append(aln)
    txt_alns.sort(key=lambda aln : aln.query_name)

    for query, group in groupby(txt_alns, lambda aln: aln.query_name):
	if query != contig:
	    continue

	alns = list(group)
	if len(alns) == 1:
	    if alns[0].is_unmapped:
		break

	    align = Alignment.from_alignedRead(alns[0], bam)
	    # find Adjacency from single gapped alignment (insertion)
	    adjs = gapped_align.find_adjs(align, contig_seq, False)
	    if adjs:
		for adj in adjs:
		    detect_itd(adj,
		               align,
		               contig_seq,
		               outdir,
		               conditions['min_len'],
		               conditions['max_apart'],
		               conditions['min_pid'],
		               debug=True)

		    # event passed by bwa mem (there is a novel insertion but failed by detect_itd()
		    # - not a duplication, will not have rna_event attribute or ITD
		    # likely a novel exon
		    if not hasattr(adj, 'rna_event') or adj.rna_event != 'ITD':
			break

		    if transcript.strand == '-':
			adj.orients = list(reversed(list(adj.orients)))

		    if debug:
			print 'itd dup', contig, adj.novel_seq, transcript.strand, adj.orients, adj.breaks

		    itds.append(adj)

	elif len(alns) == 2:
	    if alns[0].is_unmapped or alns[1].is_unmapped:
		break

	    chimeric_aligns, dubious = split_align.find_chimera(alns, 'bwa_mem', bam, min_coverage=0.8, check_haplotype=False)

	    if chimeric_aligns:
		adj = split_align.call_event(chimeric_aligns[0], chimeric_aligns[1])

		if adj.rearrangement == 'dup':
		    # make sure novel seq is None as it may be set
		    adj.novel_seq = None
		    if adj.contig_breaks[0][1] - adj.contig_breaks[0][0] > 1:
			len_changed_bases = adj.contig_breaks[0][1] - adj.contig_breaks[0][0] - 1
			adj.contig_breaks[0][0] += len_changed_bases
			adj.breaks = list(adj.breaks)
			if chimeric_aligns[0].strand == '+':
			    adj.breaks[0] += len_changed_bases
			else:
			    adj.breaks[0] -= len_changed_bases

		    elif adj.contig_breaks[0][1] - adj.contig_breaks[0][0] < 1:
			len_same_bases = adj.contig_breaks[0][0] - adj.contig_breaks[0][1] + 1
			adj.contig_breaks[0][0] -= len_same_bases
			adj.breaks = list(adj.breaks)
			if chimeric_aligns[0].strand == '+':
			    adj.breaks[0] -= len_same_bases
			else:
			    adj.breaks[0] += len_same_bases

		    if adj.breaks[0] > 0 and adj.breaks[1] > 0:
			if adj.breaks[0] < adj.breaks[1]:
			    dup_seq = txt_seq[adj.breaks[0] - 1 : adj.breaks[1]]
			else:
			    dup_seq = txt_seq[adj.breaks[1] - 1 : adj.breaks[0]]
			adj.size = len(dup_seq)
			if adj.size > 0:
			    if debug:
				print 'itd dup', contig, dup_seq, len(dup_seq)

			    if transcript.strand == '-':
				#adj.orients = list(adj.orients).reverse()
				adj.orients = list(reversed(list(adj.orients)))

			    itds.append(adj)

    if len(itds) == 1:
	exons = transcript.txt_coord_to_exon(itds[0].breaks[0]), transcript.txt_coord_to_exon(itds[0].breaks[1])
	genome_coords = transcript.txt_coord_to_genome_coord(itds[0].breaks[0]), transcript.txt_coord_to_genome_coord(itds[0].breaks[1])
	if debug:
	    print 'itd txt_to_genome', contig, exons, genome_coords, itds[0].size

	if adj_genome is not None:
	    update_adj(adj_genome,
	               transcript,
	               exons,
	               genome_coords,
	               itds[0].size,
	               itds[0].contig_breaks,
	               orients=itds[0].orients
	               )
	    itd = adj_genome
	else:
	    update_adj(itds[0],
	               transcript,
	               exons,
	               genome_coords,
	               )
	    itd = itds[0]

	return itd

    return None

def update_adj(adj, transcript, exons, genome_coords, size=None, contig_breaks=None, orients=None):
    """Updates Adjacency with all ITD-related attributes after adj is deemed to be a ITD

    Used by capture_from_transcript_align() when ITD is captured/confirmed by transcript alignment

    Args:
        adj: (Adjacency) Adjacency object
        transcript: (Transcript) Transcript object of mapped transcript
        exons: (List) exon numbers
        genome_coordinates: (List) genome coordinates of ITD
        size: (int) size of ITD
        contig_breaks: (list) contig breaks
        orients: (list) 'L','R'
    """
    adj.rna_event = 'ITD'
    adj.__class__ = Event
    adj.genes = (transcript.gene, transcript.gene)
    adj.transcripts = (transcript.id, transcript.id)
    adj.exons = exons
    adj.breaks = genome_coords
    if contig_breaks is not None:
	adj.contig_breaks = contig_breaks
    if size is not None:
	adj.size = size
    if orients is not None:
	adj.orients = orients
    adj.read_depth = 0

def detect_from_partial_alignment(align, contig_seq, transcript, ref_fasta, outdir='.', debug=False):
    """Detect ITD given a partial alignment

    1. extract the unaligned sequence
    2. see if it's duplicated within the contig sequence
    3. if it is, align against the transcript sequence to get the correct duplication coordinates

    Args:
        align: (Alignment) Alignment object
	contig_seq: (str) contig sequence - for extracting the unmapped sequence
	transcript: (Transcript) Transcript object the contig is aligned to
	ref_fasta: (FastaFile) Pysam Fasta handle to the reference genome
	outdir: (str) output directory
	debug: (Boolean) debugging mode

    Returns:
        Event object if itd is detected or None
    """
    itd = None
    entire = intspan('1-%s' % len(contig_seq))
    aligned = intspan('%s-%s' % (align.qstart, align.qend))
    for unmapped in (entire - aligned).ranges():
	unmapped_seq = contig_seq[unmapped[0] - 1 : unmapped[1]]
	dup = search_by_regex(unmapped_seq, contig_seq, 1000)

	if dup:
	    itd = confirm_with_transcript_seq(align.query, contig_seq, transcript, ref_fasta, outdir=outdir, debug=debug)

	    # transcript alignment still partial
	    if itd is None:
		if dup[1] == unmapped and dup[1][1] == len(contig_seq):
		    genome_pos = (align.qpos_to_tpos(dup[1][0] - 1),
		                  align.qpos_to_tpos(dup[0][0]))
		    itd = Event((align.target, align.target),
			        genome_pos,
			        'dup',
			        contig = align.query,
			        contig_breaks = (dup[1][0] - 1, dup[1][0]),
			        orients = ('L', 'R')
			        )
		    itd.rna_event = 'ITD'
		    itd.size = dup[1][0] - dup[0][0]
		    itd.genes = (transcript.gene, transcript.gene)
		    itd.transcripts = (transcript.id, transcript.id)
		    itd.exons = (transcript.coord_to_exon(genome_pos[0]), transcript.coord_to_exon(genome_pos[1]))

    return itd
