#!/usr/bin/env python3

###############################################################
#
#                      Singleton Rescue
#
#  Author: Nina Wang
#  Date Created: July 5, 2016
###############################################################
#  Function:
# To rescue single reads with its complimentary (SSCS/singleton) strand and enable error suppression
# - Traditionally, consensus sequences can only be made from 2 or more reads
#
# Written for Python 3.5.1
#
# Usage:
# Python3 singleton_strand_rescue.py [--singleton SingletonBAM] [--bedfile BEDFILE]
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
# 1. A BAM file containing paired singletons error corrected by its complimentary SSCS - "sscs.rescue.bam"
# 2. A BAM file containing paired singletons error corrected by its complimentary singleton - "singleton.rescue.bam"
# 3. A BAM file containing the remaining singletons that cannot be rescued as its missing a complimentary strand -
#    "rescue.remaining.bam"
# 4. A text file containing summary statistics (Total singletons, Single strand rescued singletons, % SSCS rescue,
#    Singleton strand rescued singletons, % singleton rescue, Singleton remaining (not rescued))
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

    Return consensus of complimentary reads with N for inconsistent bases.
    """
    consensus_seq = ''
    consensus_qual = []

    for i in range(read1.query_length):
        # Check to see if base at position i is the same
        if read1.query_sequence[i] == read2.query_sequence[i]:
            consensus_seq += read1.query_sequence[i]
            mol_qual = sum(read1.query_qualities[i], read2.query_qualities[i])
            # Set to max quality score if sum of qualities is greater than the threshold (Q60) imposed by genomic tools
            if mol_qual > 60:
                consensus_qual += [60]
            else:
                consensus_qual += [mol_qual]
        else:
            consensus_seq += 'N'
            consensus_qual += [0]

    return consensus_seq, consensus_qual


def strand_rescue(read_tag, duplex_tag, query_name, singleton_dict, sscs_dict=None):
    """(str, str, dict, dict) -> Pysam.AlignedSegment

    Return 'rescued' singleton read using compliment read from opposite strand (either found in SSCS or singleton).

    Quality score calculated from singleton and complimentary read. Read template retained from singleton being rescued.
    """
    read = singleton_dict[read_tag][0]

    # If SSCS bamfile provided, rescue with SSCS
    if sscs_dict is None:
        compliment_read = singleton_dict[duplex_tag][0]
    else:
        compliment_read = sscs_dict[duplex_tag][0]

    dcs = duplex_consensus(read, compliment_read)
    dcs_read = create_aligned_segment([read], dcs[0], dcs[1], query_name)

    return dcs_read


###############################
#        Main Function        #
###############################

def main():
    """Singleton rescue:
    - First rescue with SSCS bam
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
    sscs_rescue_bam = pysam.AlignmentFile('{}.sscs.rescue.bam'.format(args.singleton.split('.singleton')[0]), 'wb',
                                          template=singleton_bam)
    singleton_rescue_bam = pysam.AlignmentFile('{}.singleton.rescue.bam'.format(args.singleton.split('.singleton')[0]),
                                               'wb', template=singleton_bam)
    remaining_rescue_bam = pysam.AlignmentFile('{}.rescue.remaining.bam'.format(args.singleton.split('.singleton')[0]),
                                               'wb', template=singleton_bam)
    badRead_bam = pysam.AlignmentFile('{}.singleton.rescue.badReads.bam'.format(args.singleton.split('.singleton')[0]),
                                      "wb", template=singleton_bam)

    stats = open('{}.stats.txt'.format(args.singleton.split('.singleton')[0]), 'a')
    time_tracker = open('{}.time_tracker.txt'.format(args.singleton.split('.singleton')[0]), 'a')

    # ===== Initialize dictionaries =====
    singleton_dict = collections.OrderedDict()  # dict that remembers order of entries
    singleton_tag = collections.defaultdict(int)
    singleton_pair = collections.defaultdict(list)
    singleton_csn_pair = collections.defaultdict(list)

    sscs_dict = collections.OrderedDict()
    sscs_tag = collections.defaultdict(int)
    sscs_pair = collections.defaultdict(list)
    sscs_csn_pair = collections.defaultdict(list)

    rescue_dict = collections.OrderedDict()

    # ===== Initialize counters =====
    singleton_counter = 0
    singleton_unmapped = 0
    singleton_unmapped_mate = 0
    singleton_multiple_mappings = 0

    sscs_counter = 0
    sscs_unmapped = 0
    sscs_unmapped_mate = 0
    sscs_multiple_mappings = 0

    sscs_dup_rescue = 0
    singleton_dup_rescue = 0
    singleton_remaining = 0

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
                # Its okay to clear SSCS dicts after each region as rescue can only be done within the same coor
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
                             badRead_bam=badRead_bam,
                             duplex=None,
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
        singleton_unmapped_mate += singleton[6]
        singleton_multiple_mappings += singleton[7]

        # === Store SSCS reads in dictionaries ===
        sscs = read_bam(sscs_bam,
                        pair_dict=sscs_pair,
                        read_dict=sscs_dict,  # keeps track of paired tags
                        tag_dict=sscs_tag,
                        csn_pair_dict=sscs_csn_pair,
                        badRead_bam=badRead_bam,
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
        sscs_unmapped_mate += sscs[6]
        sscs_multiple_mappings += sscs[7]

        ######################
        #       RESCUE       #
        ######################
        for readPair in list(singleton_csn_pair.keys()):
            for tag in singleton_csn_pair[readPair]:
                counter += 1
                # Check to see if singleton can be rescued by SSCS, then by singletons
                # If not, add to 'remaining' rescue bamfile
                duplex = duplex_tag(tag)
                query_name = readPair + ':1'  # Reflect rescued singleton (non-rescues won't have our unique ID tag)

                # 1) Singleton rescue by complimentary SSCS
                if duplex in sscs_dict.keys():
                    rescue_read = strand_rescue(tag, duplex, query_name, singleton_dict, sscs_dict=sscs_dict)
                    sscs_dup_rescue += 1
                    sscs_rescue_bam.write(rescue_read)

                    del sscs_dict[duplex]
                    del singleton_dict[tag]

                # 2) Singleton rescue by complimentary Singletons
                elif duplex in singleton_dict.keys():
                    rescue_read = strand_rescue(tag, duplex, query_name, singleton_dict)
                    singleton_dup_rescue += 1
                    singleton_rescue_bam.write(rescue_read)
                    rescue_dict[tag] = duplex

                    if duplex in rescue_dict.keys():
                        del singleton_dict[tag]
                        del singleton_dict[duplex]
                        del rescue_dict[tag]
                        del rescue_dict[duplex]

                # 3) Singleton written to remaining bam if neither SSCS or Singleton duplex rescue was possible
                else:
                    remaining_rescue_bam.write(singleton_dict[tag][0])
                    singleton_remaining += 1
                    del singleton_dict[tag]

            del singleton_csn_pair[readPair]

        # Update time tracker
        time_diff = str((time.time() - start_time)/60)
        print(time_diff)
        try:
            time_tracker.write(x + ': ')
            # time_tracker.write(time_diff + '\n')
        except:
            continue

    ######################
    #       SUMMARY      #
    ######################
    sscs_rescue_frac = (sscs_dup_rescue/singleton_counter) * 100
    singleton_rescue_frac = (singleton_dup_rescue/singleton_counter) * 100

    summary_stats = '''# === Singleton Rescue ===
Total singletons: {}
SSCS strand rescued singletons: {}
% SSCS rescue: {}
Singleton strand rescued singletons: {}
% singleton rescue: {}
Singletons remaining (not rescued): {} \n'''.format(counter, sscs_dup_rescue, sscs_rescue_frac, singleton_dup_rescue, singleton_rescue_frac, singleton_remaining)

    stats.write(summary_stats)
    print(summary_stats)

    # Close files
    singleton_bam.close()
    sscs_bam.close()
    sscs_rescue_bam.close()
    singleton_rescue_bam.close()
    remaining_rescue_bam.close()
    stats.close()


###############################
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)