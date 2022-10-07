[![Release](https://img.shields.io/github/v/release/bcgsc/pavfinder?include_prereleases)](https://github.com/bcgsc/pavfinder/releases)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pavfinder/badges/latest_release_date.svg)](https://anaconda.org/bioconda/pavfinder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pavfinder/badges/installer/conda.svg)](https://anaconda.org/bioconda/pavfinder)
[![Conda](https://img.shields.io/conda/dn/bioconda/pavfinder?label=Conda)](https://anaconda.org/bioconda/pavfinder)

## *P*ost-*A*ssembly *V*ariant *Finder* (PAVFinder)

PAVFinder is a Python package that detects structural variants from *de novo* assemblies (e.g. [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom), [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss), [Trans-ABySS](http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss)).  As such, it is able to analyse both genome and transcriptome assemblies:

### genomic structural variants `pavfinder genome`
- translocations
- inversions
- duplications
- insertions
- deletions
- simple-repeat expansions/contractions

### transcriptomic structural variants `pavfinder fusion`
- gene fusions
- internal tandem duplications (ITD)
- partial tandem duplications (PTD)
- small indels
- simple-repeat expansions/contractions

### transcriptomic splice variants `pavfinder splice`
- skipped exons
- novel exons
- novel introns
- retained introns
- novel splice acceptors/donors

PAVFinder infers variants from non-contiguous (split or gapped) contig sequence alignments to the reference genome. Assemblies can be aligned to the reference genome (`c2g` alignment) using [bwa mem](http://bio-bwa.sourceforge.net/)(genome) or [gmap](http://research-pub.gene.com/gmap/)(transcriptome).  Read support for events can be gathered by aligning reads to the assembly using [bwa mem](http://bio-bwa.sourceforge.net/) (`r2c` alignment).

We provide a *T*argeted-*A*ssembly-*P*ipeline, `TAP`, to facilitate transcriptome analysis on selected genes. This requires a [multi-index Bloom Filter](http://www.bcgsc.ca/platform/bioinfo/software/biobloomtools) of targeted gene sequences to be created beforehand. Whereas whole transcriptome analysis with over 100 million read pairs can take more than 24 hours, a targeted analysis of several hundred genes (e.g. [COSMIC](http://cancer.sanger.ac.uk/census/)) can be completed within half an hour. `TAP` uses [Trans-ABySS](https://github.com/bcgsc/transabyss) for transcriptome assembly. `TAP2` is the successor of `TAP` and it uses [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom) for transcriptome assembly.

We also provide a pipeline for gene fusion detection in RNA-seq data, `Fusion-Bloom`, which couples PAVFinder with [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom). We demonstrated that it has higher senstivitiy and specificity than most state-of-the-art fusion callers.    

### Installation
See [INSTALL.md](INSTALL.md)

### Usage
See [USAGE.md](USAGE.md)

### Test Data
See [pavfinder/test](pavfinder/test) for a small dataset to test our transcriptome (`TAP`, `TAP2`, and `Fusion-Bloom`) and genome workflows.

### Publications
Readman Chiu, Ka Ming Nip, Justin Chu and Inanc Birol. **TAP: a targeted clinical genomics pipeline for detecting transcript variants using RNA-seq data**. *BMC Med Genomics* (2018) 11:79 [https://doi.org/10.1186/s12920-018-0402-6](https://doi.org/10.1186/s12920-018-0402-6)

Readman Chiu, Ka Ming Nip, Inanc Birol. **Fusion-Bloom: fusion detection in assembled transcriptomes**. *Bioinformatics* (2019) btz902 [https://doi.org/10.1093/bioinformatics/btz902](https://doi.org/10.1093/bioinformatics/btz902)
