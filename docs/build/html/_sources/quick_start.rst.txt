Quick start guide
=================
	
Installation
------------
1. Download the latest release or clone the repository::

    $ git clone https://github.com/pughlab/ConsensusCruncher.git

2. Install the dependencies: 

 ================================================ =========== =============================== 
  Program					      Version     Purpose 	                 
 ------------------------------------------------ ----------- ------------------------------- 
  `Python <https://www.python.org/>`_ 	           3.5.1   	   Run ConsensusCruncher 	     
  `BWA <http://bio-bwa.sourceforge.net/>`_ 	       0.7.15      Align reads 	  		          
  `Samtools <http://samtools.sourceforge.net/>`_   1.3.1       Sorting and indexing bamfiles  
 ================================================ =========== =============================== 

3. All required python libraries can be installed by running::

    $ pip install -r requirements.txt


Configuration
-------------
Set up ``config.ini`` with the appropriate configurations for ``fastq2bam`` and ``consensus`` modes. 
Alternatively, you can provide command-line arguments, which will overwrite ``config.ini`` parameters.

Example config file::

	[fastq2bam]
	fastq1 = # Path to FASTQ containing Read 1 of paired-end reads. [MANDATORY]
	fastq2 = # Path to FASTQ containing Read 2 of paired-end reads. [MANDATORY]
	output = # Output directory, where barcode extracted FASTQ and 
		 # BAM files will be placed in subdirectories 'fastq_tag'
		 # and 'bamfiles' respectively (dir will be created if
		 # they do not exist). [MANDATORY]
	name = # Sample name extracted from filename using read number delimiter 
		 # (e.g. Sample1_R1.fastq, delimiter = '_R', sample name = "Sample1"). [MANDATORY]
	bwa = # Path to executable BWA. [MANDATORY]
	ref = # Path to reference (BWA index). [MANDATORY]
	samtools = # Path to executable samtools. [MANDATORY]
	bpattern = # Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases). [MANDATORY]
	blist = # List of barcodes (Text file with unique barcodes on each line). [MANDATORY]
	
	# Note: You can input either a barcode list or barcode pattern or both. 
	# If both are provided, barcodes will first be matched with the list and then the constant 
	# spacer bases will be removed before the barcode is added to the header. 
	# e.g. ATNNGT means 2-bp barcode is flanked by two spacers matching 'AT' in front and 
	# 'GT' behind.

	[consensus]
	bam  = # Input BAM file with barcodes in header. [MANDATORY]
	c_output = # Output directory for consensus sequences. [MANDATORY]
	samtools = # Path to executable samtools. [MANDATORY]

	# Optional arguments
	scorrect = # Singleton correction [True or False, default: True]. 
	genome = # Genome version to determine which cytoband file to use for data splitting
		 # (hg19 or hg38, default: hg19)
	bedfile = # Bedfile, default: cytoBand.txt. WARNING: It is HIGHLY RECOMMENDED that you 
		  # use the provided cytoBand.txt unless you're working with genome build that 
		  # is not hg19 or hg38. Then a separate bedfile is needed for data segmentation. 
		  # For small BAM files, you may choose to turn off data splitting with '-b False' 
		  # and process everything all at once (Division of data is only required for 
		  # large data sets to offload the memory burden).
	cutoff = # Consensus cut-off, default: 0.7 
		 # (70% of reads must have the same base to form a consensus).
	bdelim = # Delimiter before barcode in read name 
		 # (e.g. '|' in 'HWI-D00331:196:C900FANXX:7:1110:14056:43945|TTTT')
	cleanup = # Remove intermediate files. [True or False]

UMI list example::

	TATGCGT
	AACGGAT
	AAGCCT
	AATCGCT
	ACCGAT
	ACGCAAT
	ACGTGT
	AGGACAT
	ATGTCCT
	CACAGT

	
Running ConsensusCruncher
---------------------------
ConsensusCruncher has 2 modes: 

Extract UMIs & align FASTQs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. ``fastq2bam`` extracts UMIs and aligns FASTQs with BWA to generate BAM files.
	Run ``ConsensusCruncher.py [-c CONFIG] fastq2bam`` with required input parameters 
	(see ``config.ini`` or ``--help`` for examples).
	
	This script extracts molecular barcode tags and removes spacers from unzipped FASTQ 
	files found in the input directory (file names must contain "R1" or "R2"). Barcode 
	extracted FASTQ files are written to the 'fastq_tag' directory and are subsequently 
	aligned with BWA mem. Bamfiles are written to the 'bamfile" directory under the project folder.
	
	 
Error suppression 
~~~~~~~~~~~~~~~~~
2. ``consensus`` generates single-strand and duplex consensus sequences with ‘Singleton Correction’
	Run ``ConsensusCruncher.py [-c CONFIG] consensus`` with the required input parameters 
	(see ``config.ini`` or ``--help`` for examples).
	
	Using unique molecular identifiers (UMIs), duplicate reads from the same molecule are 
	amalgamated into single-strand consensus sequences (SSCS). If complementary strands are 
	present, SSCSs can be subsequently combined to form duplex consensus sequences (DCS).
	
	If 'Singleton Correction' (SC) is enabled, single reads (singletons) can be error 
	suppressed using complementary strand. These corrected singletons can be merged with 
	SSCSs to be further collapsed into DCSs + SC.
	
	Finally, a BAM file containing only unique molecules (i.e. no duplicates) is created 
	by merging DCSs, remaining SSCSs (those that could not form DCSs), and remaining 
	singletons (those that could not be corrected).

Multiple FASTQs
---------------
``ConsensusCruncher.py`` processes one sample (2 paired-end FASTQ files or 1 BAM file) at a time. 
A sample script to generate shell scripts for multiple samples is provided `here <https://github.com/pughlab/ConsensusCruncher/blob/master/generate_scripts.sh>`_.

The script generator will create sh scripts for each file in a fastq directory.

1. The following parameters need to be changed in the config file: name, bwa, ref, samtools, 
   bpattern (alternatively if a barcode list is used instead, remove bpattern and add blist 
   as parameter). Please note: fastq1, fastq2, output, bam, and c_output can be ignored as 
   those will be updated using the generate_scripts.sh file.
2. Update generate_scripts.sh with input, output, and code_dir.
3. Run generate_scripts.sh to create sh files and then run those scripts.

Output
------
Please note the example below is for illustrative purposes only, as sample names and index 
files were removed for simplification. Order of directories and files were also altered to 
improve comprehension. ::

	.                                           Filetype
	├── sscs
	│   ├── badReads.bam                        Reads that are unmapped or have multiple alignments
	│   ├── sscs.sorted.bam                     Single-Strand Consensus Sequences (SSCS)
	│   ├── singleton.sorted.bam                Single reads (Singleton) that cannot form SSCSs
	├── sscs_SC
	|   ├── singleton.rescue.sorted.bam         Singleton correction (SC) with complementary singletons
	|   ├── sscs.rescue.sorted.bam              SC with complementary SSCSs
	|   ├── sscs.sc.sorted.bam                  SSCS combined with corrected singletons (from both rescue strategies)   [*]
	|   ├── rescue.remaining.sorted.bam         Singletons that could not be corrected
	|   ├── all.unique.sscs.sorted.bam          SSCS + SC + remaining (uncorrected) singletons
	├── dcs
	│   ├── dcs.sorted.bam                      Duplex Consensus Sequence (DCS)
	│   ├── sscs.singleton.sorted.bam           SSCSs that could not form DCSs as complementary strand was missing  
	├── dcs_SC
	│   ├── dcs.sc.sorted.bam                   DCS generated from SSCS + SC    [*]
	│   ├── sscs.sc.singleton.sorted.bam        SSCS + SC that could not form DCSs 
	│   ├── all.unique.dcs.sorted.bam           DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons
	├── read_families.txt                       Family size and frequency
	├── stats.txt                               Consensus sequence formation metrics
	├── tag_fam_size.png                        Distribution of reads across family size
	└── time_tracker.txt                        Time log

Through each stage of consensus formation, duplicate reads are collapsed together and 
single reads are written as separate files. This allows rentention of all unique molecules, 
while providing users with easy data management for cross-comparisons between error 
suppression strategies.

To simplify analyses, it would be good to focus on SSCS+SC ("sscs.sc.sorted.bam") and 
DCS+SC ("dcs.sc.sorted.bam") as highlighted above with [*].