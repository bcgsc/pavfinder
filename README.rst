PAVFinder
=========

Author: Readman Chiu (rchiu@bcgsc.ca)

Dependencies
------------

-  `BWA 0.7.4
   <http://sourceforge.net/projects/bio-bwa/files/>`_
-  `Samtools
   <http://sourceforge.net/projects/samtools/files/samtools/>`_
-  `Python 2.7.x
   <https://www.python.org/downloads/>`_
-  `Pysam-0.7.7
   <https://github.com/pysam-developers/pysam>`_
-  `intspan-0.701
   <https://pypi.python.org/pypi/intspan/>`_
-  `pybedtools-0.6.2
   <http://pythonhosted.org/pybedtools/>`_ (only if centromeric/segdup regions need to be filtered out)
-  `blastn-2.2.29
   <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/>`_ (for transcriptome analysis)

Find variants from genome assembly
----------------------------------

::

       pavfinder genome <contigs2genome.bam> bwa_mem <contig_sequence.fa> <reference_sequence.fa> <output_directory> -b <tumor_reads2contigs.bam> --normal_bam <normal_reads2contigs.bam> --min_size <minimum event size>

Find variants from transcriptome assembly
-----------------------------------------

::

       pavfinder transcriptome <contigs2genome.bam> gmap <contig_sequence.fa> <tabix_sorted_annotation.gtf> <reference_sequence.fa> <output_directory> -b <reads2contigs.bam>

