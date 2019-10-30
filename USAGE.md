* Run bash wrapper script

  ```
  pavfinder genome [`find_sv_genome.py` parameters]
  pavfinder fusion [`find_sv_transcriptome` parameters]
  pavfinder splice [`map_splice` parameters]
  ```

* Run PV to detect genomic structural variants (translocations, inversions, duplications, insertions, deletions, etc)

  ```
  find_sv_genome.py <contigs_to_genome.bam> <contigs_fasta> <genome_fasta> <outdir> --r2c <reads_to_contigs.bam>
  ```

* To generate transcripts fasta for detecting transcriptome structural variants

  ```
  extract_transcript_sequence.py <tabix-indexed gtf> <output transcripts fasta> <reference genome> --index --only_longest
  ```
  * GTF file needs to be gzipped and indexed with [Tabix](http://www.htslib.org/doc/tabix.html)
  * GTF and reference genome must have the same chromosomes
  * To create BWA index on fasta for generating contigs-to-transcripts(c2t) bam file, adding `--index`

* Run PV to detect transcriptome structural variants (fusions, read-throughs, ITDs, PTDs, InDels)

  ```
  find_sv_transcriptome.py --gbam <contigs_to_genome_bam> --tbam <contigs_to_transcripts_bam> --transcripts_fasta <indexed_transcripts_fasta> --genome_index <GMAP index genome directory and name> --r2c <reads_to_contigs_bam> <contigs_fasta> <gtf> <genome_fasta> <outdir>
  ```

* Run PV to detect novel splice variants (exon_skipping, novel_exon, novel_intron, novel_donor, novel_acceptor, retained_intron)

  ```
  map_splice.py <contigs_to_genome_bam> <contigs_fasta> <gtf> <genome_fasta> <outdir> --r2c <reads_to_contigs_bam> [--suppl_annot supplmental.gtf.gz] [--genome_bam genome.bam]
  ```

* Run full (assembly + analysis) TAP in targeted mode

  ```
  tap.py <sample> <outdir> --bf <target_genes.bf> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length>  --nprocs <number_of_processes> --params <parameters_file>
  ```

* Run full (assembly + analysis) TAP for entire transcriptome

  ```
  tap.py <sample> <outdir> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length> --nprocs <number_of_processes> --params <parameters_file>
  ```

* Run TAP for just de novo assembly

  ```
  tap.py <sample> <outdir> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length> --nprocs <number_of_processes> --only_assembly
  ```

* Run fusion-bloom for fusion calling

  ```
  source <fusion-bloom.profile>
  fusion-bloom profile=<fusion-bloom.profile> left=<fastq.gz> right=<fastq.gz> readlen=<read_length> outdir=<outdir> name=<prefix>
  ```

  example <fusion-bloom.profile>:
  ```
  export NUM_THREADS=12
  export SAMTOOLS_SORT_MEM=20G
  export GENOME=hg38
  export GMAPDB=/path/to/gmapdb_sarray/hg38
  export GTF=/path/to/gencode.v26.annotation.transcripts.sorted.gtf.gz
  export TRANSCRIPTS_FASTA=/path/to/gencode.v26.annotation.transcripts.sorted.gtf.fa
  export GENOME_FASTA=/path/to/hg38.fa
  export RNABLOOM_PARAMS='-ntcard -fpr 0.005 -chimera -extend'
  export PAVFINDER_PARAMS='--only_fusions --include_non_exon_bound_fusion --min_support 2'
  ```
