# PAVFinder

Author: Readman Chiu

## Dependencies

* [BWA 0.7.4](http://sourceforge.net/projects/bio-bwa/files/)
* [Samtools](http://sourceforge.net/projects/samtools/files/samtools/)
* [Python 2.7.x](https://www.python.org/downloads/)
* [Pysam](https://github.com/pysam-developers/pysam)
* [pybedtools](http://pythonhosted.org/pybedtools/) (only if centromeric/segdup regions need to be filtered out)

## Find structural variants

```
SV/find_sv.py <contigs2genome.bam> bwa_mem <contig_sequence.fa> <reference_sequence.fa> <output_directory> -b <tumor_reads2contigs.bam> --normal_bam <normal_reads2contigs.bam> --min_size <minimum event size>
```

More detailed documentation will come later.
