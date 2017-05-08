## Dependencies

### Software
Make sure the following software are in your PATH

- [bwa](http://bio-bwa.sourceforge.net/)
- [samtools](http://samtools.sourceforge.net/)
- [gmap](http://research-pub.gene.com/gmap/) (optional for fusion/transcriptome SV analysis)

### Python package
Install the following Python packages (e.g. using `pip insall`)

- [pysam](https://github.com/pysam-developers/pysam)
- [pybedtools](https://daler.github.io/pybedtools/)
- [intspan](https://pypi.python.org/pypi/intspan/)
- [biopython](http://biopython.org/) (optional for transcriptome analysis)

## Install PAVFinder in virtualenv

```
$ pip install virtualenv
$ virtualenv <DIR>
$ source <DIR>/bin/activate
$ pip install git+https://github.com/bcgsc/pavfinder.git#egg=pavfinder
```
