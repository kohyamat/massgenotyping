==============
massgenotyping
==============

.. image:: https://badge.fury.io/py/massgenotyping.svg
    :target: https://badge.fury.io/py/massgenotyping
    :alt: PyPI version

.. image:: https://img.shields.io/pypi/pyversions/massgenotyping.svg
    :target: https://pypi.org/project/massgenotyping
    :alt: Python versions

.. image:: https://img.shields.io/pypi/l/massgenotyping.svg
    :target: https://pypi.org/project/massgenotyping
    :alt: License


Python package for microsatellite genotyping from highly multiplexed amplicon sequencing data


Features
--------

* Semi-automatic genotyping optimized for amplicon sequencing data of microsatellite loci

* Visual genotyping with interactive plots

* Fast SSR search in sequences

* Automatic grouping and naming of alleles based on polymorphisms in both SSR and non-SSR regions

* Support for multi-core processing


Requirements
------------

* Python 3.8 or higher

* `NGmerge <https://github.com/jsh58/NGmerge>`_

* `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_

* BLASTn (included in `BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ command line applications provided by NCBI)

* Optional: `ripgrep <https://github.com/BurntSushi/ripgrep>`_


Installation
------------

This package works fine on Linux and macOS, but has not yet been fully tested on Windows.

Install from PyPI

.. code:: bash

    pip3 install massgenotyping

Install the latest version from a Git repository

.. code:: bash

    pip install git+https://github.com/kohyamat/massgenotyping


Usage
-----

.. code:: bash

    mgt [-h] SUBCOMMAND [OPTIONS]

**Subcommand list:**

* :code:`demultiplex`: demultiplex raw amplicon sequences based on primer sequences

* :code:`merge-pairs`: merge paired-end reads

* :code:`denoise`: reduce any noise that may have been generated during sequencing and PCR

* :code:`filter`: filtering for erroneous sequence variants and screening for putative alleles

* :code:`allele-check`: check allele candidates and create an allele database

* :code:`allele-call`: assign alleles to raw amplicon sequences

* :code:`show-alignment`: show a sequence alingment

The details of the options for each subcommand can be checked by :code:`mgt SUBCOMMAND -h`.


Tutorials with example data
---------------------------

Here's a step-by-step tutorial using the `example data <https://github.com/kohyamat/massgenotyping/tree/master/examples>`_.

**1. Demultiplex raw amplicon sequences based on primer sequences**

As a first step, the sequence data is split based on the primer sequence. 
The input can be one or two sequence files in the FASTQ format, or a directory containing multiple sequence files.
Primer sequences can be read from CSV or FASTA files.
Please check the example data for the format of the input data.

.. code:: bash

    mgt demultiplex examples/sequence_data -g "*_R[12]_*" -m examples/marker_data.csv

The result files are written in subdirectories within the output directory (:code:`./project` by default) for each marker.

**2. Merge paired-end reads and trim primer sequecnes**

For the paired-end sequencing data, the respective sequence pairs are merged using NGmerge program.
The following command removes the the primer sequences after merging sequence pairs.

.. code:: bash

    mgt merge-pairs ./project -m examples/marker_data.csv --trim-primer

For single-end data, this step can be skipped. The removal of the primer sequence can also be performed in the step 1.

**3. Reduce noise (optional but recommended)**

This step corrects any noise (very low-frequency point mutations) that may have been generated during sequencing or PCR.
This step is not necessarily required, but it will make the following step easier.

.. code:: bash

    mgt denoise ./project/*/*_merged.fastq.gz

**4. Filter out erroneous sequence variants**

In this step, the sequence of putative alleles is extracted for each marker in each sample,
while removing any erroneous sequence variants, such as 'stutter' sequences.
After some rough filtering, an interactive plot allows you to choose which sequence variants to keep.
You can skip this visual-checking procedure with the :code:`--force-no-visual-check` option.

.. code:: bash

    mgt filter ./project -m examples/marker_data.csv

**5. Check a multiple sequence alignment and make an allele database**

The database is created after checking the alignment of the putative allele sequences.
If necessary, you can further filter out the erroneous sequence variants.

.. code:: bash

    mgt allele-check ./project


**6. Assign alleles to raw amplicon sequences**

Finally, the following command perform a BLASTn search against the database created for each marker and assign alleles to the raw sequence data.
The genotype tables are created within the output directory.

.. code:: bash

    mgt allele-call ./project -m examples/marker_data.csv

Screenshots
-----------

.. image:: https://user-images.githubusercontent.com/6261781/78501753-205e3280-7798-11ea-98ce-32a4f631bb05.png
   :scale: 50%
   :alt: Figure 1

**Figure 1.** Checking the multiple sequence alignment across the samples (*STEP 5*).

.. image:: https://user-images.githubusercontent.com/6261781/78501825-877be700-7798-11ea-8382-3b991a42502f.png
   :scale: 50%
   :alt: Figure 2

**Figure 2.** Visual genotyping (*STEP 6*).


Contributing to massgenotyping
------------------------------

Contributions of any kind are welcome!


License
-------

`MIT <LICENSE>`_
