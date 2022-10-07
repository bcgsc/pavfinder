## Dependencies

### Software
Make sure the following software are in your `PATH`

- [bwa](http://bio-bwa.sourceforge.net/)
- [samtools](http://samtools.sourceforge.net/)
- [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [gmap](http://research-pub.gene.com/gmap/) (for Fusion-Bloom, TAP, TAP2)
- [biobloomtools](https://github.com/bcgsc/biobloom) (for TAP and TAP2)
- [RNA-Bloom](https://github.com/bcgsc/RNA-Bloom) (for Fusion-Bloom and TAP2)
- [Trans-ABySS](https://github.com/bcgsc/transabyss) (for TAP)

### Python package
Install the following Python packages (e.g. using `pip insall`)

- [pysam](https://github.com/pysam-developers/pysam)
- [pybedtools](https://daler.github.io/pybedtools/)
- [intspan](https://pypi.python.org/pypi/intspan/)
- [biopython](http://biopython.org/) (for transcriptome analysis)
- [ruffus](http://www.ruffus.org.uk/) (for TAP and TAP2)

## Install PAVFinder in a `conda` environment

```
conda create -n pavfinder
conda activate pavfinder
conda install -c bioconda pavfinder
```
It is recommended to use `mamba` instead of `conda`, i.e.
```
mamba install -c bioconda pavfinder
```

## Install PAVFinder in `virtualenv`

This is an alternative to using a `conda` environment. All required software listed above must be in your `PATH`.

```
$ pip install virtualenv
$ virtualenv <DIR>
$ source <DIR>/bin/activate
$ pip install git+https://github.com/bcgsc/pavfinder.git#egg=pavfinder
```
