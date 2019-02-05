#!/usr/bin/env python3

###############################################################
#
#      Single Stranded Consensus Sequence (SSCS) Generator
#
# Author: Nina Wang
# Date Created: Mar 24, 2016
###############################################################
# Function:
# To generate single strand consensus sequences for strand based error suppression.
# - Consensus sequence from most common base with quality score >= Q30 and greater than <cutoff> representation
# - Consensus quality score from addition of quality scores (i.e. product of error probabilities)
#
# Written for Python 3.5.1
#
# Usage:
# python3 SSCS_maker.py [--cutoff CUTOFF] [--infile INFILE] [--outfile OUTFILE] [--bedfile BEDFILE]
#
# Arguments:
# --cutoff CUTOFF     Proportion of nucleotides at a given position in a sequence required to be identical to form a
#                     consensus (Recommendation: 0.7 - based on previous literature Kennedy et al.)
#                        Example (--cutoff = 0.7):
#                           Four reads (readlength = 10) are as follows:
#                              Read 1: ACTGATACTT
#                              Read 2: ACTGAAACCT
#                              Read 3: ACTGATACCT
#                              Read 4: ACTGATACTT
#                           The resulting SSCS is: ACTGATACNT
# --infile INFILE     Input BAM file
# --outfile OUTFILE   Output BAM file
# --bedfile BEDFILE   Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt -
#                     See bed_separator.R for making your own bed file based on a target panel / specific coordinates)
#
# Inputs:
# 1. A position-sorted BAM file containing paired-end reads with duplex barcode in the header
# 2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing
#
# Outputs:
# 1. A SSCS BAM file containing paired single stranded consensus sequences - "sscs.bam"
# 2. A singleton BAM file containing single reads - "singleton.bam"
# 3. A bad read BAM file containing unpaired, unmapped, and multiple mapping reads - "badReads.bam"
# 4. A text file containing summary statistics (Total reads, Unmmaped reads, Secondary/Supplementary reads, SSCS reads,
#    and singletons) - "stats.txt"
# 5. A tag family size distribution plot (x-axis: family size, y-axis: number of reads) - "tag_fam_size.png"
# 6. A text file tracking the time to complete each genomic region (based on bed file) - "time_tracker.txt"
#
# Concepts:
#    - Read family: reads that share the same molecular barcode, genome
#                   coordinates for Read1 and Read2, cigar string, strand, flag, and read number
#    - Singleton: a read family containing only one member (a single read)
#
###############################################################

##############################
#        Load Modules        #
##############################
import pysam  # Need to install
import collections
import re
import array
from random import *
from itertools import chain
import argparse
import matplotlib.pyplot as plt
import math
import time

from consensus_helper import *


###############################
#       Helper Functions      #
###############################
def consensus_maker(readList, cutoff):
    """(list, int) -> str, list, list
    Return consensus sequence and quality score.

    Arguments:
        - readList: list of reads sharing the same unique molecular identifier
        - cutoff: Proportion of nucleotides at a given position in a sequence required to be identical to form a consensus

    Concept:
        Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position.
        - At each position, reads supporting each nucleotide is recorded along with the quality score corresponding to
          each nucleotide
        - Bases below the Phred quality cutoff (Q30) are excluded from consensus making
        - The most frequent base is added to the consensus sequence, given that the proportion of reads supporting this
          base is greater than the cutoff
        - A molecular phred quality score (consensus quality score) is determined by taking the product of errors of the
          most frequent base
        - If a majority can't be determined (i.e. a tie with 2 maximums), N will be assigned as these bases won't pass
          the proportion cut-off
    """
    # Initialize counters
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = []

    # Determine consensus for every position across read length
    readLength = readList[0].infer_query_length()

    for i in range(readLength):
        # Positions in the following lists corresponds to A, C, G, T, N
        nuc_count = [0, 0, 0, 0, 0]
        failed_nuc_count = [0, 0, 0, 0, 0]
        quality_score = [[], [], [], []]
        phred_fail = 0

        # Count bases and quality scores for position i across all reads in list
        for j in range(len(readList)):
            # Filter bases < phred quality 30 into separate list
            if readList[j].query_qualities[i] < 30:
                nuc = readList[j].query_sequence[i]
                nuc_index = nuc_lst.index(nuc)
                failed_nuc_count[nuc_index] += 1
                phred_fail += 1
            else:
                nuc = readList[j].query_sequence[i]
                nuc_index = nuc_lst.index(nuc)
                nuc_count[nuc_index] += 1
                quality_score[nuc_index].append(readList[j].query_qualities[i])

        # Find most frequent nucleotide base and quality score (don't worry about ties (2 maxes) as it won't pass the
        # proportion cut-off and N will be assigned)
        max_nuc_index = nuc_count.index(max(nuc_count))
        max_nuc_quality = quality_score[max_nuc_index]

        # Determine consensus phred quality through addition of quality scores (i.e. product of error probabilities)
        base_fail = False
        if max_nuc_quality is not []:
            mol_qual = sum(max_nuc_quality)
            # Set to max quality score if sum of qualities is greater than the threshold (Q60) imposed by genomic tools
            if mol_qual > 60:
                mol_qual = 60
        else:
            mol_qual = 0
            base_fail = True

        # Consensus only made if proportion of most common base is > cutoff (e.g. 70%)
        phred_pass_reads = len(readList) - phred_fail  # Remove number of failed bases from total count
        if phred_pass_reads != 0:
            prop_score = nuc_count[max_nuc_index]/phred_pass_reads
            if prop_score >= cutoff:
                consensus_read += nuc_lst[max_nuc_index]
                quality_consensus.append(mol_qual)
            else:
                base_fail = True
        else:
            base_fail = True

        # Set base to N if no consensus could be made
        if base_fail:
            consensus_read += 'N'
            quality_consensus.append(mol_qual)

    return consensus_read, quality_consensus


# Improve readability of argument help documentation
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


###############################
#        Main Function        #
###############################
def main():
    # Command-line parameters
    parser = ArgumentParser(formatter_class=SmartFormatter)
    parser.add_argument("--cutoff", action="store", dest="cutoff", type=float,
                        help="R|Proportion of nucleotides at a given position in a\nsequence required to be identical"
                        " to form a consensus\n(Recommendation: 0.7 - based on previous literature\nKennedy et al.)\n"
                        "   Example (--cutoff = 0.7):\n"
                        "       Four reads (readlength = 10) are as follows:\n"
                        "       Read 1: ACTGATACTT\n"
                        "       Read 2: ACTGAAACCT\n"
                        "       Read 3: ACTGATACCT\n"
                        "       Read 4: ACTGATACTT\n"
                        "   The resulting SSCS is: ACTGATACNT",
                        required=True)
    parser.add_argument("--infile", action="store", dest="infile", help="Input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="Output SSCS BAM file", required=True)
    parser.add_argument("--bdelim", action="store", dest="bdelim", default="|",
                        help="Delimiter to differentiate barcodes from read name, default: '|'")
    parser.add_argument("--bedfile", action="store", dest="bedfile",
                        help="Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt - \
                        See bed_separator.R for making your own bed file based on a target panel/specific coordinates)",
                        required=False)
    args = parser.parse_args()

    ######################
    #       SETUP        #
    ######################
    start_time = time.time()
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    # set up time tracker
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')

    # ===== Initialize dictionaries =====
    read_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    pair_dict = collections.defaultdict(list)
    csn_pair_dict = collections.defaultdict(list)

    # ===== Initialize counters =====
    unmapped = 0
    multiple_mapping = 0  # secondary/supplementary reads
    counter = 0
    singletons = 0
    SSCS_reads = 0
    bad_spacer = 0

    #######################
    #   SPLIT BY REGION   #
    #######################
    # ===== Determine data division coordinates =====
    # division by bed file if provided
    if args.bedfile is not None:
        division_coor = bed_separator(args.bedfile)
    else:
        division_coor = [1]

    # ===== Process data in chunks =====
    region = 0
    for x in division_coor:
        if division_coor == [1]:
            read_chr = None
            read_start = None
            read_end = None
        else:
            read_chr = x.split('_', 1)[0]
            read_start = division_coor[x][0]
            read_end = division_coor[x][1]

        # === Construct dictionaries for consensus making ===
        chr_data = read_bam(bamfile,
                            read_dict=read_dict,
                            tag_dict=tag_dict,
                            pair_dict=pair_dict,
                            csn_pair_dict=csn_pair_dict,
                            badRead_bam=badRead_bam,
                            duplex=None,  # this indicates bamfile is not for making DCS (thus headers are diff)
                            read_chr=read_chr,
                            read_start=read_start,
                            read_end=read_end,
                            barcode_delim=args.bdelim)

        # Set dicts and update counters
        read_dict = chr_data[0]
        tag_dict = chr_data[1]
        pair_dict = chr_data[2]
        csn_pair_dict = chr_data[3]

        counter += chr_data[4]
        unmapped += chr_data[5]
        multiple_mapping += chr_data[6]
        bad_spacer += chr_data[7]

        ######################
        #     CONSENSUS      #
        ######################
        # ===== Create consensus sequences for paired reads =====
        for readPair in list(csn_pair_dict.keys()):
            if len(csn_pair_dict[readPair]) == 2:
                for tag in csn_pair_dict[readPair]:
                    # Check for singletons
                    if tag_dict[tag] == 1:
                        singletons += 1
                        # Assign singletons our unique query name
                        read_dict[tag][0].query_name = readPair + ':' + str(tag_dict[tag])
                        singleton_bam.write(read_dict[tag][0])
                    else:
                        # Create collapsed SSCSs
                        SSCS = consensus_maker(read_dict[tag], float(args.cutoff))

                        query_name = readPair + ':' + str(tag_dict[tag])
                        SSCS_read = create_aligned_segment(read_dict[tag], SSCS[0], SSCS[1], query_name)

                        # Write consensus bam
                        SSCS_bam.write(SSCS_read)
                        SSCS_reads += 1

                    # Remove read from dictionary after writing
                    del read_dict[tag]

                # Remove key from dictionary after writing
                del csn_pair_dict[readPair]

        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')
        except:
            # When no genomic coordinates (x) provided for data division
            continue

    ######################
    #       SUMMARY      #
    ######################
    # === STATS ===
    # Note: total reads = unmapped + secondary + SSCS uncollapsed + singletons
    summary_stats = '''# === SSCS ===
Uncollapsed - Total reads: {}
Uncollapsed - Unmapped reads: {}
Uncollapsed - Secondary/Supplementary reads: {}
SSCS reads: {}
Singletons: {}
Bad spacers: {}\n'''.format(counter, unmapped, multiple_mapping, SSCS_reads, singletons, bad_spacer)

    stats.write(summary_stats)
    print(summary_stats)

    # === QC to see if there's remaining reads ===
    print('# QC: Total uncollapsed reads should be equivalent to mapped reads in bam file.')
    print('Total uncollapsed reads: {}'.format(counter))
    print('Total mapped reads in bam file: {}'.format(bamfile.mapped))

    print("QC: check dictionaries to see if there are any remaining reads")
    print('=== pair_dict remaining ===')
    if bool(pair_dict):
        for i in pair_dict:
            try:
                print(i)
                print('read remaining:')
                print(pair_dict[i][0])
                print('mate:')
                print(bamfile.mate(pair_dict[i][0]))
            except ValueError:
                print("Mate not found")
    print('=== read_dict remaining ===')
    if bool(read_dict):
        for i in read_dict:
            try:
                print(i)
                print('read remaining:')
                print(read_dict[i][0])
                print('mate:')
                print(bamfile.mate(read_dict[i][0]))
            except ValueError:
                print("Mate not found")
    print('=== csn_pair_dict remaining ===')
    if bool(csn_pair_dict):
        for i in csn_pair_dict:
            try:
                print(i)
                print(csn_pair_dict[i])
            except ValueError:
                print("Mate not found")

    # ===== write tag family size dictionary to file =====
    tags_per_fam = collections.Counter([i for i in tag_dict.values()])  # count of tags within each family size
    lst_tags_per_fam = list(tags_per_fam.items())  # convert to list [(fam, numTags)]
    with open(args.outfile.split('.sscs')[0] + '.read_families.txt', "w") as stat_file:
        stat_file.write('family_size\tfrequency\n')
        stat_file.write('\n'.join('%s\t%s' % x for x in lst_tags_per_fam))

    # ===== Create tag family size plot =====
    total_reads = sum(tag_dict.values())
    # Read fraction = family size * frequency of family / total reads
    read_fraction = [(i*j)/total_reads for i, j in lst_tags_per_fam]

    plt.bar(list(tags_per_fam), read_fraction)
    # Determine read family size range to standardize plot axis
    plt.xlim([0, math.ceil(lst_tags_per_fam[-1][0]/10) * 10])
    plt.savefig(args.outfile.split('.sscs')[0]+'_tag_fam_size.png')

    # ===== Close files =====
    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    badRead_bam.close()


###############################
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)

