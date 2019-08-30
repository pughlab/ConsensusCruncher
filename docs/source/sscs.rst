Single-strand consensus sequences (SSCS)
========================================

``SSCS_maker.py``

**Function:**
To generate single strand consensus sequences for strand based error suppression.

	- Consensus sequence from most common base with quality score >= Q30 and greater than <cutoff> representation
	- Consensus quality score from addition of quality scores (i.e. product of error probabilities)

(Written for Python 3.5.1)

**Usage:**
	python3 SSCS_maker.py [--cutoff CUTOFF] [--infile INFILE] [--outfile OUTFILE] 
	[--bedfile BEDFILE]

**Arguments:**

+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --cutoff CUTOFF    | Proportion of nucleotides at a given position in a sequence required to be identical to form a consensus                                                                                 |
|                    | 	- Recommendation: 0.7 based on previous literature Kennedy et al.                                                                                                                       |
|                    | 	- Example (--cutoff = 0.7) - four reads (readlength = 10) are as follows:                                                                                                               |
|                    | 		- Read 1: ACTGATACTT                                                                                                                                                            |
|                    | 		- Read 2: ACTGAAACCT                                                                                                                                                            |
|                    | 		- Read 3: ACTGATACCT                                                                                                                                                            |
|                    | 		- Read 4: ACTGATACTT                                                                                                                                                            |
|                    | 	- The resulting SSCS is: ACTGATACNT                                                                                                                                                     |
+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --infile INFILE    | Input BAM file                                                                                                                                                                           |
+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --outfile OUTFILE  | Output BAM file                                                                                                                                                                          |
+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| --bedfile BEDFILE  | Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt)													|
+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


**Inputs:**
	1. A position-sorted BAM file containing paired-end reads with duplex barcode in the header
	2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing

**Outputs:**
	1. A SSCS BAM file containing paired single stranded consensus sequences - "sscs.bam"
	2. A singleton BAM file containing single reads - "singleton.bam"
	3. A bad read BAM file containing unpaired, unmapped, and multiple mapping reads - "badReads.bam"
	4. A text file containing summary statistics (Total reads, Unmmaped reads, Secondary/Supplementary reads, SSCS reads,
	   and singletons) - "stats.txt"
	5. A tag family size distribution plot (x-axis: family size, y-axis: number of reads) - "tag_fam_size.png"
	6. A text file tracking the time to complete each genomic region (based on bed file) - "time_tracker.txt"

**Concepts:**
   - Read family: reads that share the same molecular barcode, genome coordinates for 
     Read1 and Read2, cigar string, strand, flag, and read number
   - Singleton: a read family containing only one member (a single read)
