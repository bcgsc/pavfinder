v0.1.1

- fixed up variant size-range filtering so that both `min_size` and `max_size` can be used at the same time
- will get rid of any inversion with size = 1
- read_support now tallied by `Variant` instead of `Adjacency` so that an inversion with 2 breakpoints/adjacencies will have its read_support (of each adjacency) summed up (bwa_mem, no multi-mapping) for filtering

v0.1.2

- Fixed bug in ignoring orientations when creating fusion
- new functions in generating BED track by `alignment.py`:
  - output correct strand (useful for displaying strand-specific contig alignment)
  - minimum size `--min_size` can be specified so that small contig alignments can be skipped
- Fixed bug in mis-labelling reciprocal translocations as insertions
- Don't use Pysam AlignedRead.rlen for checking if sequence is chimeric as BWA versions later than 0.7.4 outputs only chimeric portion of sequence other than entire sequence

v0.2.0 (all transcriptome changes)

- reports coverage/depth of all exon-exon junctions assembled in BED format `junctions.bed`
- reports coverage/depth of reference 5'`ref5_jn_depth` and 3'`ref5_jn_depth` junctions in `events.tsv` for gene fusions and all novel splicing events
- renamed header "spanning_reads" to "support_reads"
- accepts supplementary gene annotation GTF `--suppl_annot` for checking novelty of splicing events
- will not call ITD on homopolymer expansion
- changed event-label of same-gene chimera to actual rearrangement (e.g. a duplication within the same gene will be called 'dup' but not 'fusion')

v0.3.0

- added support for detecting potential transposable element insertion and repeat expansion and contraction
- changed support read detection workflow
- added test data

v0.4.0

- integrated [pavfinder_transcriptome](https://github.com/bcgsc/pavfinder_transcriptome)
- changed defaults of `subseq_len` and `probe_len` in `find_sv_transcriptome.py` for detecting events in test data
- changed `adjacencies.tsv` output from genome to `adjacencies.bedpe` to agree with transcriptome output
- changed event label `repeat_reduction` to `repeat_contraction`
- increased maximum microhomology from 5 to 10 for transcriptome analysis

v0.4.1

- chimera subseq checking improvements/fixes: determine subseq length from fasta file instead of cigar string in BAM file (hard-clip creates problem), multi-map checking requires mapping of entire sequence instead of just checking targets

v0.4.2

- minor changes to `tap.py`:
  - requires `samtools` version >= 1.0 because `samtools sort` command changes
  - creates ruffus history file in output directory
  - explicity index fasta file after merging is done and before analysis scripts are run
  - determine `num_proc` parameter for running `find_sv.py`
  - changed name of one of duplicated `format_read_pairs()` functions
  - fixed bug if suppl_annot is not str

v0.4.3

- changes to handle novel untemplated sequence at ITD breakpoint
	- uses `blastn` to check for duplication, and doesn't require end-to-end matching, but requires insertion sequence to be at least 15bp(hard-coded) long; allows bases at edges un-matched (specified by `max_novel_length` default=10)
	- when parsing partial alignments, allow clipped bases (specified by `max_novel_length` again)
	- introduce parameter `min_dup_size`: minimum insertion size to check if it's a duplication (default: 15)
	- changed default `--max_novel_len` from 20 to 10
- added `--only_fusions` parameter to only look for fusions
- fixed bug in detecting retained_intron: retained_intron in last match block wasn't detected.

v0.4.4

- bugfix: if retained_intron happened in first and only alignment block, the variable 'j' will be exposed as un-initialized. Initialize it to 0
- always assume every alignment will have blocks but there are actually cases where Pysam cannot determine alignment blocks. Capture and skip those cases.
- fixed r2c bug in tap when running in whole transcriptome mode
- updated expected output files for tap
- updated refgenes.sorted.gtf.gz and setup.py to copy ensGenes.sorted.gtf.gz

v1.0.0
- added `fusion-bloom`(Make script) for running [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom) instead of Trans-ABySS
- modified `tap.py` to deal with multimatch headers from BBT 2.1.1
- modified `tap.py` to parse in assembly parameters from cfg file
- bugfix: extract gene names from c2t.bam even when no adjacency/event is identified

v1.1.0
- allows 1 base to differ from canonical splice site (GTAG) in detecting novel donors or acceptors because of possibility of splice-site mutation
- can input `genome_bam` to detect support level of splice site variant
- search through main annotation to make sure novel event is novel; before only supplementary annotations is used for this purpose, and when no supplementary annotation is provided, an erroneous "novel" event may be reported because wrong transcript is assigned for contig

v1.2.0
- associate novel splicing events with genomic variants given WGS alignments
- bugfix: retained_intron captured by multiple contigs not reported because `link` attribute is lost after merging adjacencies - fixed
- bugfix: successive retained introns not captured - fixed
- bugfix: remove retained_intron called if encompassing exon (not necessarily having exact boundaries) is found in annotation
- remove redundant novel_acceptor and novel_donor calls of novel_introns

v1.3.0
- added `tap2.py`, using [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom) instead of Trans-ABySS for assembly

v1.4.0
- allow fusion flanking_pair support reads to overlap span window (default size 8bp, 4 on either side of breakpoint) so that flanking_pair support of short fragment can also be captured
- changed r2c to run BWA-mem in paired-end mode in `tap2.py`
- fixed bug in `tap2.py` extracting `RNA-Bloom` version when there is stderr message while running `rnabloom -v'

v1.5.0
- added `--min_overhang` (default=4) for changing minimum overhang size for gathering spanning reads from r2c in `find_sv_transcriptome.py`
