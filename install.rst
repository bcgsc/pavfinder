Installation
------------

Install the following softwares and put there executables in ``PATH``:

-  `BWA 0.7.4
   <http://sourceforge.net/projects/bio-bwa/files/>`_
-  `Samtools
   <http://sourceforge.net/projects/samtools/files/samtools/>`_
-  `GMAP 2014-01-21
   <http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-01-21.tar.gz>`_
-  `blastn-2.2.29
   <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/>`_ 

Install PAVFinder:

::

   $ pip install virtualenv
   $ virtualenv <DIR>
   $ source <DIR>/bin/activate
   $ pip install git+https://github.com/bcgsc/pavfinder.git#egg=pavfinder
