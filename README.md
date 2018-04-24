# README - Tools for Duplex Sequencing #

This README would normally document whatever steps are necessary to get your application up and running.

## What is this repository for? ##

* Runs the molecular barcode / duplex sequencing pipeline for duplicate removal and error suppression

## How do I get set up? ##

### Configuration ###
Set up DuplexPipeline.sh
- Setup cluster configuration (default SGE cluster with highmem.q) 
Note: depending on the size of your bamfiles, software might require a lot of memory resources

~~~~
sh [Duplex Sequencing git directory]/DuplexPipeline.sh Input_dir Output_dir 
~~~~
Although the DuplexPipeline.sh script takes in a bedfile, this needs to be a specially formatted bedfile (using the bed\_separator.R tool). 
It is HIGHLY recommended you run the script without a bedfile (if this is your time), so the default "cytoband.txt" will be used to separate the bam file for processing. 

Scripts will be created for each bamfile found in the input directory. 
Each bamfile will be processed through the following:
1. Single stranded consensus sequence (SSCS)
2. Duplex consensus sequence (DCS)
3. Singleton Correction (SC)
4. Merge SSCS and Singleton Correction (SSCS_SC)	
5. DCS from SSCS + Singleton Correction (DCS_SC)
6. Generate all unique molecular bamfiles 

### Dependencies ###
This pipeline requires the following dependencies:

| Program | Version | Purpose                                    |
| ------- | ------- | ------------------------------------------ |
| Python3 | 3.5.1   | Consensus sequence pipeline                |
| Numpy   | 1.11.0  | Python library for scientific computing    |
| Pandas  | 0.19.2  | Python library for data analysis           |
| Pysam   | 0.9.0   | Python interface for working with bamfiles |
| Samtools| 1.3.1   | Sorting and indexing bamfiles              |
| Picard  | 2.6.0   | Merging bamfiles                           |
| Java    | 8       | Used with Picard to merge bamfiles         |

# Intro to Molecular Barcoding #
![Scheme](beta/script_overview.png =500x500)

**Duplex sequencing schematic:** 
An uncollapsed BAM file is first processed through SSCS_maker.py to create an error suppressed single-strand 
consensus sequences (SSCS) BAM file and an uncorrected Singleton BAM file. The single reads can be corrected
through singleton_correction.py, which rescues singletons with its complementary SSCS or singleton. SSCS 
reads can be directly made into duplex consensus sequences (DCS) or merged with corrected singletons to create
an expanded pool of DCS reads (Figure illustrates singleton correction merged work flow).

### Bamfiles ###
* Uncollapsed: Original bamfiles
* SSCS: Single strand consensus sequences
* SSCS_SC: Single strand consensus sequences with Singleton Correction
* SSCS_SC_Singletons: SSCS + SC that could not be made into DCSs
* All_unique_sscs: Single strand consensus sequences + Singleton Correction + remaining (unrescued) singletons 
* DCS: Duplex consensus sequences
* DCS_SC: Duplex consensus sequences from SSCS_SC
* All_unique_dcs: Duplex consensus sequences from SSCS_SC + SSCS_SC_Singletons + remaining singletons
* Singletons: Single reads

### Who do I talk to? ###

* Repo owner or admin (Nina)
