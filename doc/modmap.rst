===========================================================
modmap: software for analyzing modified bases in the genome
===========================================================
:Homepage: http://github.com/hesselberthlab/modmap
:Author: Jay R. Hesselberth
:Organization: University of Colorado School of Medicine
:Address: Department of Biochemistry and Molecular Genetics,
          Aurora CO 80045 USA
:Copyright: 2013
:Last updated: |today|
:Version: |version|

.. currentmodule:: modmap

The modmap software is a simple set of analysis tools supporting methods
for mapping nucleobases in genomes.

Programs
========

- genome_nuc_freqs: calculates nucleotide frequencies from a genome FASTA
  file and reports frequencies in a tabular format.

- origin_analysys: calculates signals relative to annotations of
  replication origins.

- random_dist: assesses distribution of data in genome. outputs tables of
  observed and randomized data for downstream comparison

- summary_table: generates table of observed counts at positions across
  samples. uses bedtools.multi_intersect.


