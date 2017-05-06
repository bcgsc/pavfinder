import os
from setuptools import setup, find_packages
from pavfinder import __version__

setup(
    name='pavfinder',
    version=__version__,
    description='Post Assembly Variant Finder',
    long_description='Identifies genomic structural variants or transcriptomic splice variants given the alignments of assembly contigs',
    url='https://github.com/bcgsc/pavfinder.git',
    author='Readman Chiu',
    author_email='rchiu@bcgsc.ca',
    license='BCCA',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(),
    install_requires = [
        'pysam>=0.8.1',
        'pybedtools>=0.6.9',
        'intspan>=0.701',
        'numpy>=1.9.2',
        'pandas',
        'biopython',
        ],
    package_data = {'pavfinder': ['test/*', 'scripts/*']},
    data_files = [('config', ['pavfinder/cfg/tap.cfg']), ('test', ["pavfinder/test/*"])],
    scripts = ['pavfinder/scripts/pavfinder',
               'pavfinder/scripts/check_support.py',
               'pavfinder/scripts/find_sv_genome.py',
               'pavfinder/scripts/find_sv_transcriptome.py',
               'pavfinder/scripts/map_splice.py',
               'pavfinder/scripts/extract_transcript_sequence.py',
               'pavfinder/scripts/tap.py',
               'pavfinder/scripts/rescue_fusion.py',
               ],
)
