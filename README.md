# ConsensusCruncher #

ConsensusCruncher is a tool that suppresses errors in next-generation sequencing data by using unique molecular identifers (UMIs) to amalgamate reads derived from the same DNA template into a consensus sequence.

## Quick start ##
### Dependencies ###
This pipeline requires the following dependencies:

| Program | Version | Purpose                                    |
| ------- | ------- | ------------------------------------------ |
| [Python3](https://www.python.org/) | 3.5.1   | Consensus sequence pipeline                |
| [Numpy](http://www.numpy.org/)   | 1.11.0  | Python library for scientific computing    |
| [Pandas](https://pandas.pydata.org/)  | 0.19.2  | Python library for data analysis           |
| [Pysam](https://pypi.org/project/pysam/#description)   | 0.9.0   | Python interface for working with bamfiles |
| [Samtools](http://samtools.sourceforge.net/)| 1.3.1   | Sorting and indexing bamfiles              |
| [Picard](https://broadinstitute.github.io/picard/)  | 2.6.0   | Merging bamfiles                           |
| Java    | 8       | Used with Picard to merge bamfiles         |

### Configuration ###
Set up fastq_to_bam.sh and ConsensusCruncher.sh with appropriate configurations:
1. **Cluster**: Current settings are set-up for Sun Grid Engine HPC clusters. Depending 
on the size of your bamfiles, high memory nodes may be required. 
2. **Modules**: Update the scripts if you are using a different version of any of the 
required programs (Python, Samtools, Picard, Java).

### Running ConsensusCruncher ###
1. Run **fastq_to_bam.sh** with required input parameters: \
-i  Input directory [MANDATORY] \
-o  Output project directory [MANDATORY] \
-p  Project name [MANDATORY] \
-r  Reference (BWA index) [MANDATORY] \
-b  Barcode length [MANDATORY] \
-s  Spacer length [MANDATORY] 
    
This script extracts molecular barcode tags and removes spacers from unzipped FASTQ
files found in the input directory (file names must contain "R1" or "R2"). Barcode
extracted FASTQ files are written to the 'fastq_tag' directory and are subsequently
aligned with BWA mem. Bamfiles are written to the 'bamfile" directory under the
project folder.

2. Run **ConsensusCruncher.sh** with the required input parameters: \
-i  Input bamfile directory [MANDATORY] \
-o  Output project directory [MANDATORY]
  
This script amalgamates duplicate reads in bamfiles into single-strand consensus
sequences (SSCS), which are subsequently combined into duplex consensus sequences
(DCS). Singletons (reads lacking duplicate sequences) are corrected, combined
with SSCS to form SSCS + SC, and further collapsed to form DCS + SC. Finally,
files containing all unique molecules (a.k.a. no duplicates) are created for SSCS
and DCS.

## How it works ##
Unique molecular identifiers (UMIs) composed of molecular barcodes and sequence features are used aggregate reads derived from the same strand of a template molecule. Amalgamation of such reads into single strand consensus sequences (SSCS) removes discordant bases, which effectively eliminates polymerase and sequencer errors. Complementary SSCSs can be subsequently combined to form a duplex consensus sequence (DCS), which eliminates asymmetric strand artefacts such as those that develop from oxidative damage. 

Conventional UMI-based strategies rely on redundant sequencing from both template strands to form consensus sequences and cannot error suppress single reads (singleton). We enable singleton correction using complementary duplex reads in the absence of redundant sequencing. 

## Overview ##
<img src="https://user-images.githubusercontent.com/13406244/39268149-03b4c12a-489d-11e8-8011-f85ec8a82f39.png" width="50%" height="50%">

**ConsensusCruncher schematic:**
* An uncollapsed bamfile is first processed through SSCS_maker.py to create an error-suppressed single-strand 
consensus sequence (SSCS) bamfile and an uncorrected singleton bamfile. 
* The singletons can be corrected through singleton_correction.py, which error suppress singletons with its complementary SSCS or singleton read. 
* SSCS reads can be directly made into duplex consensus sequences (DCS) or merged with corrected singletons to create
an expanded pool of DCS reads (Figure illustrates singleton correction merged work flow).


## Example ##
In order to create consensus sequences, we first need to process fastq files into bam files. Sample fastq files can be found under the test folder.

### Fastqs to Bams ###
Given **fastq** as the input directory, *fastq_to_bam.sh* removes the spacer region and extracts the barcode tag from each sequencing read into the header with *extract_barcode.py*.

```
PATH="[insert path to ConsensusCruncher repo]"
BWAPATH="[insert path to BWA index]"
sh fastq_to_bam.sh -i ./$PATH/ConsensusCruncher/test/fastq -o ./$PATH/ConsensusCruncher/test -b 2 -s 1 -f T -r ./$BWAPATH/BWAIndex/genome.fa
```

In the sample dataset, we utilized 2-bp barcodes and 1-bp spacers. While the barcodes for each read can be one of 16 possible combinations (4^2), the spacer is an invariant "T" base used to ligate barcodes onto each end of a DNA fragment. Thus, a spacer filter (-f) should be imposed to remove faulty reads. Barcodes from read 1 and read 2 are extracted and combined together before being added to the header. 

```
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
@HWI-D00331:196:C900FANXX:5:1101:1332:2193|ATGG/1
AAGCCCCAGGCAGTTGCTAATGATGGGAGCTTAGTGCACAAGGGCTGGGCCTCCCTCTTGGAGCTGAACATTGTTTCTTGGGGACGGCTGTGCCCACCTCAGCGGGGAGGCAAGGATTAAATC
+
CCGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGEGGGGGBGGGGGGGGGGGGGGGGGGGGGGGEGG1:FGFGGGGGGGGG/CB>DG@GGGGGGG<DGGGGAAGGEGGB>DGGGEGGG/@G

Read2:
@HWI-D00331:196:C900FANXX:5:1101:1332:2193|ATGG/2
GGGCTCCAGCCCTGATTTCCTCCCCCAGCCCTGCAGGGCTCAGGTCCAGAGGACACAAGTTTAACTTGCGGGTGGTCACTTGCCTCGTGCGGTGACGCCATGGTGCCCTCTCTGTGCAGCGCA
+
BCGGGGEGGGGFGGGGGGGGGGGGGGGGGGGGGGB:FCGGGGGGGGGGEGGGGGGGG=FCGG:@GGGEGBGGGAGFGDE@FGGGGGFGFGEGDGGGFCGGDEBGGGGGGGEG=EGGGEEGGG#

```

FASTQ files with extracted barcodes are placed in the **fastq_tag** directory and are subsequently aligned with BWA to generate BAMs in the **bamfiles** folder.

```
. 
├── bamfiles 
├── fastq
├── fastq_tag
└── qsub
```

### ConsensusCruncher ###
*ConsensusCruncher.sh* creates a **consensus** directory containing folders for each bam file. Separate bash scripts are generated for for each bamfile using the workflow illustrated above. 

Bamfiles are generated and grouped according to type of error suppression (SSCS vs DCS) and whether Singleton Correction (SC) was implemented. 
```
. 
├── bamfiles 
├── consensus 
│   ├── LargeMid_56_L005 
│   │   ├── dcs 
│   │   ├── dcs_SC 
│   │   ├── sscs 
│   │   └── sscs_SC 
... 
│   ├── LargeMid_62_L006
│   │   ├── dcs
│   │   ├── dcs_SC
│   │   ├── sscs
│   │   └── sscs_SC
│   └── qsub
├── fastq
├── fastq_tag
└── qsub
```
Within a sample directory (e.g. LargeMid_56_L005), you will find the following files:

Please note the example below is for illustrative purposes only, as sample names and index files were removed for simplification. Order of directories and files were also altered to improve comprehension.
```
.                                           Filetype
├── sscs
│   ├── badReads.bam                        Reads that are unmapped or have multiple alignments
│   ├── sscs.sorted.bam                     Single-Strand Consensus Sequences (SSCS)
│   ├── singleton.sorted.bam                Single reads (Singleton) that cannot form SSCSs
├── sscs_SC
|   ├── singleton.rescue.sorted.bam         Singleton correction (SC) with complementary singletons
|   ├── sscs.rescue.sorted.bam              SC with complementary SSCSs
|   ├── sscs.sc.sorted.bam                  SSCS combined with corrected singletons (from both rescue strategies)   [*]
|   ├── rescue.remaining.sorted.bam         Singletons that could not be corrected
|   ├── all.unique.sscs.sorted.bam          SSCS + SC + remaining (uncorrected) singletons
├── dcs
│   ├── dcs.sorted.bam                      Duplex Consensus Sequence (DCS)
│   ├── sscs.singleton.sorted.bam           SSCSs that could not form DCSs as complementary strand was missing  
├── dcs_SC
│   ├── dcs.sc.sorted.bam                   DCS generated from SSCS + SC    [*]
│   ├── sscs.sc.singleton.sorted.bam        SSCS + SC that could not form DCSs 
│   ├── all.unique.dcs.sorted.bam           DCS (from SSCS + SC) + SSCS_SC_Singletons + remaining singletons
├── read_families.txt                       Family size and frequency
├── stats.txt                               Consensus sequence formation metrics
├── tag_fam_size.png                        Distribution of reads across family size
└── time_tracker.txt                        Time log

```
Through each stage of consensus formation, duplicate reads are collapsed together and single reads are written as separate files. This allows rentention of all unique molecules, while providing users with easy data management for cross-comparisons between error suppression strategies. 

To simplify analyses, it would be good to focus on SSCS+SC ("sscs.sc.sorted.bam") and DCS+SC ("dcs.sc.sorted.bam") as highlighted above with [*].


### Who do I talk to? ###
* Repo owner or admin (Nina T. Wang)
