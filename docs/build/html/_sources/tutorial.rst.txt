Tutorial 
========

Sample FASTQ files can be found under the test folder. Please note these FASTQ are only 
for testing purposes. For the full FASTQs used in our paper, please download the data from 
the NCBI Sequence Read Archive (SRA; https://www.ncbi.nlm.nih.gov/sra/) under access 
numbers SRP140497 and SRP141184.

In order to create consensus sequences, we first need to process FASTQ files into BAM files. 

FASTQs to BAMs
--------------
Given FASTQs as input files, ``fastq2bam`` mode removes the spacer region and extracts the 
barcode tag from each sequencing read into the header with extract_barcode.py. The tag
removed FASTQs are then aligned with BWA mem into BAM files (Arguments can be provided in 
the `config.ini` file or as command-line arguments. Please note command-line arguments 
over-writes config.ini arguments). ::

	REPO="[insert path to ConsensusCruncher repo]"
	BWAPATH="[insert path to BWA]"
	BWAINDEX="[insert path to BWA INDEX]"
	BWAPATH="[insert path to SAMTOOLS]"

	python ConsensusCruncher.py fastq2bam --fastq1 $REPO/test/fastq/LargeMid_56_L005_R1.fastq 
	--FASTQ2 $REPO/test/fastq/LargeMid_56_L005_R2.fastq -o $REPO/test -b $BWAPATH -r $BWAIndex 
	-s $SAMTOOLS -bpattern NNT 

In the sample dataset, we utilized 2-bp (NN) barcodes and 1-bp (T) spacers. While the 
barcodes for each read can be one of 16 possible combinations (4^2), the spacer is an 
invariant "T" base used to ligate barcodes onto each end of a DNA fragment. Thus, a spacer 
filter is imposed to remove faulty reads. Barcodes from read 1 and read 2 are extracted and 
combined together before being added to the header. ::

	READ FROM SEQUENCER
	Read1:
	@HWI-D00331:196:C900FANXX:5:1101:1332:2193 1:N:0:ACGTCACA   [<-- HEADER]
	ATTAAGCCCCAGGCAGTTGCTAATGATGGGAGCTTAGTGCACAAGGGCTGGGCCTCCCTCTTGGAGCTGAACATTGTTTCTTGGGGACGGCTGTGCCCACCTCAGCGGGGAGGCAAGGATTAAATC  [<-- SEQUENCE]
	+
	BCCCCGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGEGGGGGBGGGGGGGGGGGGGGGGGGGGGGGEGG1:FGFGGGGGGGGG/CB>DG@GGGGGGG<DGGGGAAGGEGGB>DGGGEGGG/@G  [<-- QUALITY SCORE]

	Read2:
	@HWI-D00331:196:C900FANXX:5:1101:1332:2193 2:N:0:ACGTCACA
	GGTGGGCTCCAGCCCTGATTTCCTCCCCCAGCCCTGCAGGGCTCAGGTCCAGAGGACACAAGTTTAACTTGCGGGTGGTCACTTGCCTCGTGCGGTGACGCCATGGTGCCCTCTCTGTGCAGCGCA
	+
	BBBBCGGGGEGGGGFGGGGGGGGGGGGGGGGGGGGGGB:FCGGGGGGGGGGEGGGGGGGG=FCGG:@GGGEGBGGGAGFGDE@FGGGGGFGFGEGDGGGFCGGDEBGGGGGGGEG=EGGGEEGGG#

	------

	AFTER BARCODE EXTRACTION AND SPACER ("T") REMOVAL
	Read1:
	@HWI-D00331:196:C900FANXX:5:1101:1332:2193|AT.GG/1
	AAGCCCCAGGCAGTTGCTAATGATGGGAGCTTAGTGCACAAGGGCTGGGCCTCCCTCTTGGAGCTGAACATTGTTTCTTGGGGACGGCTGTGCCCACCTCAGCGGGGAGGCAAGGATTAAATC
	+
	CCGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGEGGGGGBGGGGGGGGGGGGGGGGGGGGGGGEGG1:FGFGGGGGGGGG/CB>DG@GGGGGGG<DGGGGAAGGEGGB>DGGGEGGG/@G

	Read2:
	@HWI-D00331:196:C900FANXX:5:1101:1332:2193|AT.GG/2
	GGGCTCCAGCCCTGATTTCCTCCCCCAGCCCTGCAGGGCTCAGGTCCAGAGGACACAAGTTTAACTTGCGGGTGGTCACTTGCCTCGTGCGGTGACGCCATGGTGCCCTCTCTGTGCAGCGCA
	+
	BCGGGGEGGGGFGGGGGGGGGGGGGGGGGGGGGGB:FCGGGGGGGGGGEGGGGGGGG=FCGG:@GGGEGBGGGAGFGDE@FGGGGGFGFGEGDGGGFCGGDEBGGGGGGGEG=EGGGEEGGG#

FASTQ files with extracted barcodes are placed in the fastq_tag directory and are 
subsequently aligned with BWA to generate BAMs in the bamfiles folder. ::

	. 
	├── bamfiles 
	├── fastq
	├── fastq_tag
	└── qsub

ConsensusCruncher
-----------------
``consensus`` mode creates a consensus directory and folders for each bam file. BAM files 
undergo consensus construction through the workflow illustrated above. Output BAMs are 
grouped according to type of error suppression (SSCS vs DCS) and whether Singleton 
Correction (SC) was implemented. ::

	. 
	├── bamfiles 
	├── consensus 
	│   ├── LargeMid_56_L005 
	│   │   ├── dcs 
	│   │   ├── dcs_SC 
	│   │   ├── sscs 
	│   │   └── sscs_SC 
	... 
	│   ├── LargeMid_62_L006
	│   │   ├── dcs
	│   │   ├── dcs_SC
	│   │   ├── sscs
	│   │   └── sscs_SC
	│   └── qsub
	├── fastq
	├── fastq_tag
	└── qsub

Within a sample directory (e.g. LargeMid_56_L005), you will find the following files:

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

Within the stats file you should expect to see the following (Please note as this is a 
test dataset, the number of consensus reads is very low)::

	# === SSCS ===
	Uncollapsed - Total reads: 19563
	Uncollapsed - Unmapped reads: 17
	Uncollapsed - Secondary/Supplementary reads: 24
	SSCS reads: 0
	Singletons: 19522
	Bad spacers: 0

	# QC: Total uncollapsed reads should be equivalent to mapped reads in bam file.
	Total uncollapsed reads: 19563
	Total mapped reads in bam file: 19563
	QC: check dictionaries to see if there are any remaining reads
	=== pair_dict remaining ===
	=== read_dict remaining ===
	=== csn_pair_dict remaining ===
	0.02919737100601196
	# === DCS ===
	SSCS - Total reads: 0
	SSCS - Unmapped reads: 0
	SSCS - Secondary/Supplementary reads: 0
	DCS reads: 0
	SSCS singletons: 0 

	# === Singleton Correction ===
	Total singletons: 19522
	Singleton Correction by SSCS: 0
	% Singleton Correction by SSCS: 0.0
	Singleton Correction by Singletons: 4
	% Singleton Correction by Singletons : 0.020489703923778302
	Uncorrected Singletons: 19518 

	0.020557292302449546
	# === DCS - Singleton Correction ===
	SSCS SC - Total reads: 4
	SSCS SC - Unmapped reads: 0
	SSCS SC - Secondary/Supplementary reads: 0
	DCS SC reads: 2
	SSCS SC singletons: 0 