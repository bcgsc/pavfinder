v0.1.1

- fixed up variant size-range filtering so that both ``min_size`` and ``max_size`` can be used at the same time
- will get rid of any inversion with size = 1
- read_support now tallied by ``Variant`` instead of ``Adjacency`` so that an inversion with 2 breakpoints/adjacencies will have its read_support (of each adjacency) summed up (bwa_mem, no multi-mapping) for filtering

v0.1.2

- Fixed bug in ignoring orientations when creating fusion
- new functions in generating BED track by ``alignment.py``:
  
  + output correct strand (useful for displaying strand-specific contig alignment)
  + minimum size (``--min_size``) can be specified so that small contig alignments can be skipped

- Fixed bug in mis-labelling reciprocal translocations as insertions
- Don't use Pysam AlignedRead.rlen for checking if sequence is chimeric as BWA versions later than 0.7.4 outputs only chimeric portion of sequence other than entire sequence
