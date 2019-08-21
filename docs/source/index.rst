ConsensusCruncher - Extended Error Suppression
==============================================
Welcome to the ConsensusCruncher documentation! 

ConsensusCruncher is a tool that suppresses errors in next-generation sequencing data by 
using unique molecular identifiers (UMIs) to amalgamate reads derived from the same DNA 
template into a consensus sequence.

To learn more about ConsensusCruncher and its applications, see our  publication in
`Nucleic Acids Research <https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz474/5498633>`_. 

.. toctree::
	:maxdepth: 1
	:caption: Contents:
	
	install
	quick_start
# 	tutorial

Modes
-----
ConsensusCruncher has 2 modes: 

* ``fastq2bam`` extracts UMIs and aligns FASTQs with BWA to generate BAM files.
	This script extracts molecular barcode tags and removes spacers from unzipped FASTQ 
	files found in the input directory (file names must contain "R1" or "R2"). Barcode 
	extracted FASTQ files are written to the 'fastq_tag' directory and are subsequently 
	aligned with BWA mem. Bamfiles are written to the 'bamfile" directory under the 
	project folder.
* ``consensus`` generates single-strand and duplex consensus sequences with 'Singleton Correction'
	This script amalgamates duplicate reads in bamfiles into single-strand consensus 
	sequences (SSCS), which are subsequently combined into duplex consensus sequences 
	(DCS). Singletons (reads lacking duplicate sequences) are corrected, combined with 
	SSCS to form SSCS + SC, and further collapsed to form DCS + SC. Finally, files 
	containing all unique molecules (a.k.a. no duplicates) are created for SSCS and DCS.
	
.. toctree::
   :maxdepth: 1
   :caption: Commands:
   SSCS_maker
   Singleton_correction
   DCS_maker
