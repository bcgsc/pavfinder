```
python find_sv.py test.bam bwa_mem test.fa /path/of/bwa_indexed/hg19.fa /path/of/output_directory --min_size 10
```
`test.fa` contains sequences of:
* 2 insertion events of which the 2 breakpoints are captured in separate sequences
* 1 insertion event of which the 2 breakpoints are captured in the same sequence
* 1 reciprocal transclation event captured by 2 contig sequences
* 1 tri-nucleotide repeat-expansion event
