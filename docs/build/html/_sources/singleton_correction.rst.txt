Singleton Correction
====================

``singleton_correction.py``

**Function:**
To correct single reads with its complementary (SSCS/singleton) strand and enable error suppression
	- Traditionally, consensus sequences can only be made from 2 or more reads

(Written for Python 3.5.1)

**Usage:**
Python3 singleton_correction.py [--singleton Singleton BAM] [--bedfile BEDFILE]

**Arguments:**

+--------------------------+-----------------------------------------------------------------------------------------+
| --singleton SingletonBAM | input singleton BAM file                                                                |
+--------------------------+-----------------------------------------------------------------------------------------+
| --bedfile BEDFILE        | Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt) |
+--------------------------+-----------------------------------------------------------------------------------------+

**Inputs:**
	1. A position-sorted BAM file containing paired-end single reads with barcode identifiers in the header/query name
	2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing

**Outputs:**
	1. A BAM file containing paired singletons error corrected by its complementary SSCS - "sscs.correction.bam"
	2. A BAM file containing paired singletons error corrected by its complementary singleton - "singleton.correction.bam"
	3. A BAM file containing the remaining singletons that cannot be corrected as its missing a complementary strand -
	   "uncorrected.bam"
	4. A text file containing summary statistics (Total singletons, Singleton Correction by SSCS, % Singleton Correction by SSCS,
	   Singleton Correction by Singletons, % Singleton Correction by Singletons, Uncorrected Singletons)
	   - "stats.txt" (Stats pended to same stats file as SSCS)

**Concepts:**
   - Read family: reads that share the same molecular barcode, chr, and start coordinates for Read1 and Read2
   - Singleton: single read with no PCR duplicates (family size = 1)
