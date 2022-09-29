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
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [biopython](http://biopython.org/) (for transcriptome analysis)
- [ruffus](http://www.ruffus.org.uk/) (for TAP)

## Set up a `conda` environment

It is recommended to use `mamba` instead of `conda`.

```
conda create -n pavfinder
conda activate pavfinder
conda install --file conda_requirements.txt -c conda-forge -c bioconda
```
```
mamba create -n pavfinder
mamba activate pavfinder
mamba install --file conda_requirements.txt -c conda-forge -c bioconda
```

Before running `tap`, `tap2`, `fusion-bloom`, or other scripts, set your `PATH` and `PYTHONPATH` after activating your `conda` environment, e.g.
```
conda activate pavfinder
export PATH=/path/to/pavfinder/scripts:${PATH}
export PYTHONPATH=/path/to/pavfinder
```

## Install PAVFinder in `virtualenv`

This is an alternative to using a `conda` environment. All required software listed above must be in your `PATH`.

```
$ pip install virtualenv
$ virtualenv <DIR>
$ source <DIR>/bin/activate
$ pip install git+https://github.com/bcgsc/pavfinder.git#egg=pavfinder
```
