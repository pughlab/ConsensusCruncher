Duplex consensus sequences (DCS)
================================

``DCS_maker.py``

**Function:**
To generate duplex/double-strand consensus sequences for molecule based error suppression.
	- Consensus sequence from bases with Phred quality >= 30
	- Consensus quality score from addition of quality scores (i.e. product of error probabilities)

(Written for Python 3.5.1)

**Usage:**
Python3 DCS_maker.py [--infile INFILE] [--outfile OUTFILE] [--bedfile BEDFILE]

**Arguments:**

+-------------------+-----------------------------------------------------------------------------------------+
| --infile INFILE   | input BAM file                                                                          |
+-------------------+-----------------------------------------------------------------------------------------+
| --outfile OUTFILE | output BAM file                                                                         |
+-------------------+-----------------------------------------------------------------------------------------+
| --bedfile BEDFILE | Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt) |
+-------------------+-----------------------------------------------------------------------------------------+

**Inputs:**
	1. A position-sorted BAM file containing paired-end reads with SSCS consensus identifier in the header/query name
	2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing

**Outputs:**
	1. A BAM file containing paired double stranded consensus sequences - "dcs.bam"
	2. A SSCS singleton BAM file containing SSCSs without reads from the complementary strand - "sscs.singleton.bam"
	3. A text file containing summary statistics (Total SSCS reads, Unmmaped SSCS reads, Secondary/Supplementary SSCS
	   reads, DCS reads, and SSCS singletons) - "stats.txt" (Stats pended to same stats file as SSCS)

**Concepts:**
   - Read family: reads that share the same molecular barcode, chr, and start coordinates for Read1 and Read2
   - SSCS Singleton: a SSCS read without its complementary read