Data is provided for testing the detection of genome structural variants, transcriptome structural and splice variants, and the TAP targeted pipeline.  

## pavfinder genome
* input `test/genome`
  * `test.fa`: 
    * 2 sequences corresponding to a reciprocal translocation event, 
    * 5-kb insertions where only breakpoints are captured,
    * 1 250bp insertion where the entire insertion sequence is captured,
    * 1 3.6kb deletion
    * 1 13kb duplication
    * 2 sequences corresponding to a 860bp inversion
    * 1 tri-nucleotide length polymorphism
  * `test.bam`: bwa mem alignment of `test.fa` to hg19

* output `test/genome/expected_output`
  * `adjacencies.bedpe`, `variants.vcf`

```
pavfinder genome test.bam test.fa /path/to/hg19.fa /path/to/output_directory --min_size 10
```

## pavfinder transcriptome
* input `test/transcriptome/`
  * `test.fa`:
    * one NUP98/NSD1 fusion
    * one FLT3 internal tandem duplication (ITD)
    * one KMT2A partial tandem duplication (PTD)
    * one NPM1 insertion
  * `c2g.bam`: GMAP alignment of `test.fa` to hg19
  * `c2t.bam`: BWA mem alignment of `test.fa` to `refGene.fa`
  * `r2c.bam`: BWA mem alignment of reads to `test.fa`
  * `refGene.fa`: refseq transcript sequences which `test.fa` align to
  * `refGene.sorted.gtf.gz`: refseq GTF file sorted and indexed by tabix

* output`test/transcriptome/expected_output/pavfinder`:
  * `sv.bedpe`: fusions, ITD, PTD, indel events
  * `novel_splicing.bed`: no novel splice events detected in given data
  * `junctions.bed`: all exon junctions with read-depth
  * `mappings.tsv`: gene/transcript coverage by each contig sequence

```
pavfinder fusion --gbam c2g.bam --tbam c2t.bam --transcripts_fasta refGene.fa --genome_index /path/to/gmapdb hg19 --r2c r2c.bam test.fa refGene.sorted.gtf.gz /path/to/hg19.fa /path/to/output_directory
pavfinder splice c2g.bam test.fa refGene.sorted.gtf.gz /path/to/hg19.fa /path/to/output_directory --r2c r2c.bam
```

## tap
* input `test/transcriptome`
  * `test_1.fastq.gz`, `test_2.fastq.gz`
  * `test.cfg`: 
    * full paths to `genome_fasta` and `genome_index`(GMAP hg19 index) must be specified
    * full paths to downloaded `refGene.fa` and `refGene.sorted.gtf.gz` can be specified for `transcripts_fasta` and `gtf`
  
* output `test/transcriptome/expected_output/tap.tar.gz`

```
python tap.py test /path/to/output_directory --bf test_genes.bf --fq test_1.fastq.gz test_2.fastq.gz --k 32 62 92 --readlen 100 --params test.cfg --remove_fq
```
