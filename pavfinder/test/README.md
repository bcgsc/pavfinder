## pavfinder genome
```
python find_sv_genome.py test.bam bwa_mem test.fa /path/to/hg19.fa /path/to/output_directory --min_size 10
```

## pavfinder transcriptome
```
python find_sv_transcriptome.py 
--gbam c2g.bam --tbam c2t.bam --transcripts_fasta refGene.fa --genome_index /path/to/gmapdb hg19 --r2c r2c.bam test.fa refGene.sorted.gtf.gz /path/to/hg19.fa /path/to/output_directory
```

## tap
```
python tap.py test /path/to/output_directory --bf test_genes.bf --fq test_1.fastq.gz test_2.fastq.gz --k 32 62 92 --readlen 100 --params test.cfg --remove_fq
```
