Extract barcodes
================

``extract_barcodes.py``

**Function:**
To isolate duplex barcodes from paired-end sequence reads and store in FASTQ headers after removal of spacer regions.

(Written for Python 3.5.1)

**USAGE:** 
	python3 extract_barcodes.py [--read1 READ1] [--read2 READ2] [--outfile OUTFILE] 
	[--blen BARCODELEN] [--slen SPACERLEN] [--sfilt SPACERFILT]

Arguments:

+---------------------+--------------------------------------------------------------------------+
| --read1 READ1       | Input FASTQ file for Read 1 (unzipped)                                   |
+---------------------+--------------------------------------------------------------------------+
| --read2 READ2       | Input FASTQ file for Read 2 (unzipped)                                   |
+---------------------+--------------------------------------------------------------------------+
| --outfile OUTFILE   | Output FASTQ files for Read 1 and Read 2 using given filename            |
+---------------------+--------------------------------------------------------------------------+
| --bpattern BPATTERN | Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases) |
+---------------------+--------------------------------------------------------------------------+
| --blist BARCODELIST | List of correct barcodes                                                 |
+---------------------+--------------------------------------------------------------------------+

**Barcode design:**
	N = random / barcode bases
	
	A | C | G | T = Fixed spacer bases
	
	e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front and 'GT' behind

**Inputs:**
	1. A FASTQ file containing first-in-pair (Read 1) reads
	2. A FASTQ file containing second-in-pair (Read 2) reads

**Outputs:**
	1. A Read 1 FASTQ file with barcodes added to the FASTQ header
	2. A Read 2 FASTQ file with barcodes added to the FASTQ header
	3. A text file summarizing barcode stats