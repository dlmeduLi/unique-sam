Unique-Sam
==========

**Unique-Sam** is a simple command line tool to remove the duplicated
alignments in the `SAM <https://github.com/samtools/hts-specs>`__ file.
If the MAPQ field of the alignment is available, *unique-sam* will keep
one and only one alignment with the highest score. Otherwise,
*unique-sam* will calculate a score according to the alignment's MD or
CIGAR field and use the calculated value to remove the duplicated
alignments.

Install
=======

-  Install with the source code, in the source folder:

   .. code:: bash

       python setup.py install

-  If you have `**pip** <https://pip.pypa.io/en/latest/index.html>`__
   installed, you can simply run

   .. code:: bash

       pip install unique-sam

   After installation you can access **unique-sam** from your command
   line.

Usage
=====

**unique-sam** need a SAM format file to run properly. In your command
line environment:

.. code:: bash

    unique-sam input.sam -o output.sam

For more about **unique-sam** run:

.. code:: bash

    unique-sam --help

Unique Strategy
===============

Following strategies are applied to find the unique & the best alignment

1. Keep the alignment pair that has the highest score. If more than one
   pairs are found to have the same "Highest Score", these pairs will be
   removed.
2. Read1 and Read2 should be mapped on different strands.
3. The segment length decided by the read pairs should be longer than
   0.7 \* read length

Log File
========

All removed alignments will be written into log file ``input.sam.log``
under current folder. Each line of the log file start with a symbol and
followed by the deleted alignment (the original alignment record in the
``input.sam``). The symbol describe the reason of why this/these
alignments should be removed. The specification of these symbols are
listed in the follow table:

\| Symbol \| Description \| \| ------ \| ----------- \| \| **!** \|
Error lines \| \| **<** \| Low score alignments \| \| **=** \| Pairs
with more than one best score \| \| **~** \| Read pair mapped on the
same strand \| \| **?** \| Segment length too short \| \| --- \|
---------------------- \|

Copyright
=========

Copyright (c) 2015 dlmeduLi@163.com
