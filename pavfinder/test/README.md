Data is provided for testing the detection of genome structural variants, transcriptome structural and splice variants, TAP, and Fusion-Bloom pipelines.  

## pavfinder genome
* input `genome`
  * `test.fa`: 
    * 2 sequences corresponding to a reciprocal translocation event, 
    * 2 5-kb insertions where only breakpoints are captured,
    * 1 250bp insertion where the entire insertion sequence is captured,
    * 1 3.6kb deletion
    * 1 13kb duplication
    * 2 sequences corresponding to a 860bp inversion
    * 1 tri-nucleotide length polymorphism
  * `c2g.bam`: bwa mem alignment of `test.fa` to hg19
  * `r2c.bam`: bwa mem alignment of reads to `test.fa`

* output `genome/expected_output`
  * `adjacencies.bedpe`, `variants.vcf`: events detected from contig sequences
  * `adjacencies_filtered.bedpe`, `variants_filtered.vcf`: events with read support
  * `coords.tsv`: intermediate file listing event contig coordinates for gathering read support
  * `support.tsv`: intermediate results file of read support

```
pavfinder genome c2g.bam test.fa /path/to/hg19.fa /path/to/output_directory --min_size 10 --r2c r2c.bam
```

## pavfinder transcriptome
* input `transcriptome`
  * `test.fa`:
    * 1 NUP98/NSD1 fusion
    * 1 EIF4E3-FOXP1 read_through
    * 1 FLT3 internal tandem duplication (ITD)
    * 1 KMT2A partial tandem duplication (PTD)
    * 1 NPM1 insertion
    * 3 CEBPA deletions
    * 1 ATXN3 trinucleotide-repeat expansion polymorphism
    * 1 CEBPA trinucleotide-repeat contraction polymorphism
    * various FLT3 novel splicing events (skipped exons, novel splice donors/acceptors)
  * `c2g.bam`: GMAP alignment of `test.fa` to hg19
  * `c2t.bam`: BWA mem alignment of `test.fa` to `refGene.fa`
  * `r2c.bam`: BWA mem alignment of reads to `test.fa`
  * `refGene.fa`: refseq transcript sequences which `test.fa` align to
  * `refGene.sorted.gtf.gz`: refseq GTF file sorted and indexed by tabix
  * `acembly.sorted.gtf.gz`: aceview GTF file as supplementary annotation for novel splicing detection

* output`transcriptome/expected_output/pavfinder`:
  * `sv.bedpe`: fusions, ITD, PTD, indel events
  * `novel_splicing.bed`: no novel splice events detected in given data
  * `junctions.bed`: all exon junctions with read-depth
  * `mappings.tsv`: gene/transcript coverage by each contig sequence

```
pavfinder fusion --gbam c2g.bam --tbam c2t.bam --transcripts_fasta refGene.fa --genome_index /path/to/gmapdb hg19 --r2c r2c.bam test.fa refGene.sorted.gtf.gz /path/to/hg19.fa /path/to/output_directory
pavfinder splice c2g.bam test.fa refGene.sorted.gtf.gz /path/to/hg19.fa /path/to/output_directory --r2c r2c.bam
```

## tap and tap2
* `tap2` uses RNA-Bloom as the transcriptome assembler instead of Trans-ABySS, and runs the same way as `tap` (make sure rnabloom is in $PATH)
* input `transcriptome`
  * `test_1.fastq.gz`, `test_2.fastq.gz`
  * `test.cfg`: 
    * specify full paths to `genome_fasta` and `genome_index`(GMAP hg19 index)
    * specify full paths to `transcripts_fasta`(`refGene.fa` provided) and `gtf`(`refGene.sorted.gtf.gz` provided)
  * `test_genes.bf`: Bloom filter corresponding to the genes described above containing the various events
* output `transcriptome/expected_output/tap.tar.gz` or `transcriptome/expected_output/tap2.tar.gz`

The Bloom filter provided for testing with `tap` and `tap2` was generated with the following spaced seeds:
```
biobloommimaker -F -b 0.8 -t 64 -p test_genes -S "001100111101110000001100101001111001110110011010010001101001 111011000000001001100001011110000010011011111101111110010010 010010011111101111110110010000011110100001100100000000110111 100101100010010110011011100111100101001100000011101111001100" ./test_genes/*.fa
```

```
tap test /path/to/output_directory --bf test_genes.bf --fq test_1.fastq.gz test_2.fastq.gz --k 32 62 92 --readlen 100 --params test.cfg --remove_fq
tap2 test /path/to/output_directory --bf test_genes.bf --fq test_1.fastq.gz test_2.fastq.gz --readlen 100 --params test.cfg --remove_fq
```

## fusion-bloom
* input `transcriptome`
  * `test_1.fastq.gz`, `test_2.fastq.gz`
  * `test.profile`
	* specify full paths to `GMAPDB` and `GENOME_FASTA`(GMAP hg19 index and hg19 fasta file)
	* specify full paths to `TRANSCRIPTS_FASTA`(`refGene.fa` provided) and `GTF`(`refGene.sorted.gtf.gz` provided)
* output `transcriptome/expected_output/fusion-bloom.tar.gz`

```
fusion-bloom profile=test.profile left=test_1.fastq.gz right=test_2.fastq.gz readlen=100 name=test
```
