#!/usr/bin/env python3

###############################################################
#
#                      Singleton Correction
#
#  Author: Nina Wang
#  Date Created: July 5, 2016
###############################################################
# Function:
# To correct single reads with its complementary (SSCS/singleton) strand and enable error suppression
# - Traditionally, consensus sequences can only be made from 2 or more reads
#
# Written for Python 3.5.1
#
# Usage:
# Python3 singleton_correction.py [--singleton Singleton BAM] [--bedfile BEDFILE]
#
# Arguments:
# --singleton SingletonBAM  input singleton BAM file
# --bedfile BEDFILE         Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt -
#                           See bed_separator.R for making your own bed file based on specific coordinates)
#
# Inputs:
# 1. A position-sorted BAM file containing paired-end single reads with barcode identifiers in the header/query name
# 2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing
#
# Outputs:
# 1. A BAM file containing paired singletons error corrected by its complementary SSCS - "sscs.correction.bam"
# 2. A BAM file containing paired singletons error corrected by its complementary singleton - "singleton.correction.bam"
# 3. A BAM file containing the remaining singletons that cannot be corrected as its missing a complementary strand -
#    "uncorrected.bam"
# 4. A text file containing summary statistics (Total singletons, Singleton Correction by SSCS, % Singleton Correction by SSCS,
#    Singleton Correction by Singletons, % Singleton Correction by Singletons, Uncorrected Singletons)
#    - "stats.txt" (Stats pended to same stats file as SSCS)
#
# Concepts:
#    - Read family: reads that share the same molecular barcode, chr, and start
#                   coordinates for Read1 and Read2
#    - Singleton: single read with no PCR duplicates (family size = 1)
#
###############################################################

##############################
#        Load Modules        #
##############################
import pysam  # Need to install
import collections
from argparse import ArgumentParser
import time
import math
import os
import inspect

from consensus_helper import *


###############################
#       Helper Functions      #
###############################
def duplex_consensus(read1, read2):
    """(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) ->
    pysam.calignedsegment.AlignedSegment
    Return consensus of complementary reads with N for inconsistent bases.
    """
    consensus_seq = ''
    consensus_qual = []

    for i in range(read1.query_length):
        # Check to see if base at position i is the same
        if read1.query_sequence[i] == read2.query_sequence[i] and \
                        read1.query_qualities[i] > 29 and read2.query_qualities[i] > 29:
            consensus_seq += read1.query_sequence[i]
            mol_qual = sum([read1.query_qualities[i], read2.query_qualities[i]])
            # Set to max quality score if sum of qualities is greater than the threshold (Q60) imposed by genomic tools
            if mol_qual > 60:
                consensus_qual += [60]
            else:
                consensus_qual += [mol_qual]
        else:
            consensus_seq += 'N'
            consensus_qual += [0]

    return consensus_seq, consensus_qual


def strand_correction(read_tag, duplex_tag, query_name, singleton_dict, sscs_dict=None):
    """(str, str, dict, dict) -> Pysam.AlignedSegment
    Return 'corrected' singleton using complement read from opposite strand (either found in SSCS or singleton).

    Quality score calculated from singleton and complementary read. Read template based on singleton.
    """
    read = singleton_dict[read_tag][0]

    # If SSCS bamfile provided, corrected with SSCS
    if sscs_dict is None:
        complement_read = singleton_dict[duplex_tag][0]
    else:
        complement_read = sscs_dict[duplex_tag][0]

    dcs = duplex_consensus(read, complement_read)
    dcs_read = create_aligned_segment([read], dcs[0], dcs[1], query_name)

    return dcs_read


###############################
#        Main Function        #
###############################

def main():
    """Singleton correction:
    - First correct with SSCS bam
    - Rescue remaining singletons with singleton bam
    """
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--singleton", action="store", dest="singleton", help="input singleton BAM file",
                        required=True, type=str)
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
    singleton_bam = pysam.AlignmentFile(args.singleton, "rb")
    # Infer SSCS bam from singleton bamfile (by removing extensions)
    sscs_bam = pysam.AlignmentFile('{}.sscs{}'.format(args.singleton.split('.singleton')[0],
                                                      args.singleton.split('.singleton')[1]), "rb")
    sscs_correction_bam = pysam.AlignmentFile('{}.sscs.correction.bam'.format(args.singleton.split('.singleton')[0]), 'wb',
                                          template=singleton_bam)
    singleton_correction_bam = pysam.AlignmentFile('{}.singleton.correction.bam'.format(args.singleton.split('.singleton')[0]),
                                               'wb', template=singleton_bam)
    uncorrected_bam = pysam.AlignmentFile('{}.uncorrected.bam'.format(args.singleton.split('.singleton')[0]),
                                               'wb', template=singleton_bam)

    stats = open('{}.stats.txt'.format(args.singleton.split('.singleton')[0]), 'a')

    # ===== Initialize dictionaries =====
    singleton_dict = collections.OrderedDict()  # dict that remembers order of entries
    singleton_tag = collections.defaultdict(int)
    singleton_pair = collections.defaultdict(list)
    singleton_csn_pair = collections.defaultdict(list)

    sscs_dict = collections.OrderedDict()
    sscs_tag = collections.defaultdict(int)
    sscs_pair = collections.defaultdict(list)
    sscs_csn_pair = collections.defaultdict(list)

    correction_dict = collections.OrderedDict()

    # ===== Initialize counters =====
    singleton_counter = 0
    singleton_unmapped = 0
    singleton_multiple_mappings = 0

    sscs_counter = 0
    sscs_unmapped = 0
    sscs_multiple_mappings = 0

    sscs_dup_correction = 0
    singleton_dup_correction = 0
    uncorrected_singleton = 0

    counter = 0  # Total singletons

    #######################
    #   SPLIT BY REGION   #
    #######################
    if args.bedfile is not None:
        division_coor = bed_separator(args.bedfile)
    else:
        division_coor = [1]

    last_chr = "chrM"
    for x in division_coor:
        if division_coor == [1]:
            read_chr = None
            read_start = None
            read_end = None
        else:
            read_chr = x.split('_', 1)[0]
            read_start = division_coor[x][0]
            read_end = division_coor[x][1]

            # === Reset dictionaries ===
            if last_chr != read_chr:
                singleton_tag = collections.defaultdict(int)
                # Its okay to clear SSCS dicts after each region as correction can only be done within the same coor
                sscs_dict = collections.OrderedDict()
                sscs_tag = collections.defaultdict(int)
                sscs_pair = collections.defaultdict(list)
                sscs_csn_pair = collections.defaultdict(list)

                last_chr = read_chr

        # === Store singleton reads in dictionaries ===
        singleton = read_bam(singleton_bam,
                             pair_dict=singleton_pair,
                             read_dict=singleton_dict,  # keeps track of paired tags
                             tag_dict=singleton_tag,
                             csn_pair_dict=singleton_csn_pair,
                             badRead_bam=None,
                             duplex=True,
                             read_chr=read_chr,
                             read_start=read_start,
                             read_end=read_end
                             )

        singleton_dict = singleton[0]
        singleton_tag = singleton[1]
        singleton_pair = singleton[2]
        singleton_csn_pair = singleton[3]

        singleton_counter += singleton[4]
        singleton_unmapped += singleton[5]
        singleton_multiple_mappings += singleton[6]

        # === Store SSCS reads in dictionaries ===
        sscs = read_bam(sscs_bam,
                        pair_dict=sscs_pair,
                        read_dict=sscs_dict,  # keeps track of paired tags
                        tag_dict=sscs_tag,
                        csn_pair_dict=sscs_csn_pair,
                        badRead_bam=None,
                        duplex=True,
                        read_chr=read_chr,
                        read_start=read_start,
                        read_end=read_end
                        )

        sscs_dict = sscs[0]
        sscs_tag = sscs[1]
        sscs_pair = sscs[2]
        sscs_csn_pair = sscs[3]

        sscs_counter += sscs[4]
        sscs_unmapped += sscs[5]
        sscs_multiple_mappings += sscs[6]

        ########################
        # Singleton Correction #
        ########################
        for readPair in list(singleton_csn_pair.keys()):
            for tag in singleton_csn_pair[readPair]:
                counter += 1
                # Check to see if singleton can be corrected by SSCS, then by singletons
                # If not, add to 'uncorrected' bamfile
                duplex = duplex_tag(tag)
                query_name = readPair + ':1'  # Reflect corrected singleton (uncorrected won't have our unique ID tag)

                # 1) Singleton correction by complementary SSCS
                if duplex in sscs_dict.keys():
                    corrected_read = strand_correction(tag, duplex, query_name, singleton_dict, sscs_dict=sscs_dict)
                    sscs_dup_correction += 1
                    sscs_correction_bam.write(corrected_read)

                    del sscs_dict[duplex]
                    del singleton_dict[tag]

                # 2) Singleton correction by complementary Singletons
                elif duplex in singleton_dict.keys():
                    corrected_read = strand_correction(tag, duplex, query_name, singleton_dict)
                    singleton_dup_correction += 1
                    singleton_correction_bam.write(corrected_read)
                    correction_dict[tag] = duplex

                    if duplex in correction_dict.keys():
                        del singleton_dict[tag]
                        del singleton_dict[duplex]
                        del correction_dict[tag]
                        del correction_dict[duplex]

                # 3) Singleton written to remaining bam if neither SSCS or Singleton duplex correction was possible
                else:
                    uncorrected_bam.write(singleton_dict[tag][0])
                    uncorrected_singleton += 1
                    del singleton_dict[tag]

            del singleton_csn_pair[readPair]


    ######################
    #       SUMMARY      #
    ######################
    sscs_correction_frac = (sscs_dup_correction/singleton_counter) * 100
    singleton_correction_frac = (singleton_dup_correction/singleton_counter) * 100

    summary_stats = '''# === Singleton Correction ===
Total singletons: {}
Singleton Correction by SSCS: {}
% Singleton Correction by SSCS: {}
Singleton Correction by Singletons: {}
% Singleton Correction by Singletons : {}
Uncorrected Singletons: {} \n'''.format(counter, sscs_dup_correction, sscs_correction_frac, singleton_dup_correction, singleton_correction_frac, uncorrected_singleton)

    stats.write(summary_stats)
    print(summary_stats)

    # Close files
    singleton_bam.close()
    sscs_bam.close()
    sscs_correction_bam.close()
    singleton_correction_bam.close()
    uncorrected_bam.close()
    stats.close()


###############################
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)