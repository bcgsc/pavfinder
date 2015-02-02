import os
from setuptools import setup, find_packages
execfile(os.path.dirname(os.path.realpath(__file__)) + "/scripts/version.py")
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
        'pysam>=0.7.7,<0.8.1',
        'pybedtools>=0.6.2',
        'intspan>=0.701',
        ],
	scripts = ['pavfinder/scripts/pavfinder',
		],
)
