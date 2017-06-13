#!/usr/bin/env python3

###############################################################
#
#                          DCS Maker
#
# Author: Nina Wang
# Last Modified: May 18, 2016
# Date Created: Mar 24, 2016
###############################################################
# Function:
# Written for Python 3.5.1
#
# Inputs:
# 1. A position-sorted paired-end BAM file containing reads with a duplex tag
#    in the header.
#
# Outputs:
# 1. A paired-end BAM file containing single stranded consensus sequences.
# 2. A text file containing summary statistics (DCS reads)
#    - Read family: reads that share the same molecular barcode, chr, and start
#                   coordinates for Read1 and Read2
#    - Singleton: a read family containing only one member (a single read)
#
# usage: DCS_maker.py [--infile INFILE]
#                      [--outfile OUTFILE]
# optional arguments:
# --infile INFILE     input BAM file
# --outfile OUTFILE   output BAM file
#
# To run script: python3 $cwd/DCS_maker.py --infile $identifier.sscs.bam --outfile $identifier.dcs.bam
#
# python3 DCS_maker.py --infile MEM-001_KRAS.sscs.bam --outfile MEM-001_KRAS.dcs.bam
#
###############################################################

import pysam # Need to install
import collections
import re
import array
from random import randint
from argparse import ArgumentParser
import math

from consensus_helper import *

###############################
#       Helper Functions      #
###############################


def dcs_consensus_tag(tag, ds):
    '''(str, str) -> str
    Return consensus tag for duplex reads.

    >>> dcs_consensus_tag('TTCA_7_55259315_7_55259454_98M_98M_neg:3', 'CATT_7_55259315_7_55259454_98M_98M_pos:6')
    'CATT_TTCA_7_55259315_7_55259454_98M_98M'

    >>> dcs_consensus_tag('CTTC_23_74804535_23_74804611_98M_98M_neg:2', 'TCCT_23_74804535_23_74804611_98M_98M_pos:3')
    'CTTC_TCCT_23_74804535_23_74804611_98M_98M'

    >>> dcs_consensus_tag('TTTC_7_140477735_7_140477790_98M_98M_neg:3', 'TCTT_7_140477735_7_140477790_98M_98M_pos:2')
    'TCTT_TTTC_7_140477735_7_140477790_98M_98M'
    '''
    barcode = tag.split('_')[0]
    duplex_barcode = ds.split('_')[0]
    tag_coor = tag.split('_', 1)[1].rsplit('_', 1)[0]

    dcs_query_name = "{}_{}_{}".format(min(barcode, duplex_barcode),
                                       max(barcode, duplex_barcode),
                                       tag_coor)

    return dcs_query_name


def duplex_consensus(read1, read2):
    '''(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> pysam.calignedsegment.AlignedSegment

    Return consensus of 2 reads with N for variant bases.
    '''
    consensus_seq = ''
    qual_consensus = []

    for i in range(read1.query_length):
        if read1.query_sequence[i] == read2.query_sequence[i]:
            consensus_seq += read1.query_sequence[i]
            P = 10**(-(read1.query_qualities[i]/10)) * 10**(-(read2.query_qualities[i]/10))
            Q = round(-10 * math.log10(P))
            if Q > 62:
                qual_consensus += [62]
            else:
                qual_consensus += [Q]
        else:
            consensus_seq += 'N'
            qual_consensus += [0]

    return consensus_seq, qual_consensus


###############################
#        Main Function        #
###############################

def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True, type = str)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output BAM file", required = True, type = str)
    parser.add_argument("--bedfile", action = "store", dest="bedfile", help="output BAM file", required = False, type = str)
    args = parser.parse_args()

    start_time = time.time()

    args.infile = str(args.infile)
    args.outfile = str(args.outfile)

    SSCS_bam = pysam.AlignmentFile(args.infile, "rb")
    DCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = SSCS_bam)
    SSCS_singleton = pysam.AlignmentFile('{}.sscs.singleton.bam'.format(args.outfile.split('.dcs')[0]), "wb", template = SSCS_bam)
    badRead_bam = pysam.AlignmentFile('{}.dcs.badReads.bam'.format(args.outfile.split('.dcs')[0]), "wb", template = SSCS_bam)

    stats = open('{}.stats.txt'.format(args.outfile.split('.dcs')[0]), 'a')
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.dcs')[0]), 'a')

    # ===== Initialize dictionaries and convert SSCS data into dict=====
    read_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    pair_dict = collections.defaultdict(list)
    csn_pair_dict = collections.defaultdict(list)

    unmapped = 0
    unmapped_mate = 0
    multiple_mapping = 0  # secondary/supplementary reads
    counter = 0
    SSCS_singletons = 0
    multiple_mappings = 0

    duplex_count = 0
    duplex_dict = collections.defaultdict(int)

    # ===== Determine data division coordinates =====
    # division by bed file if provided
    if args.bedfile is not None:
        division_coor = bed_separator(args.bedfile)
    else:
        division_coor = [1]

    # ===== Process data in chunks =====
    for x in division_coor:
        if division_coor == [1]:
            read_chr = None
            read_start = None
            read_end = None
        else:
            read_chr = x.rsplit('_', 1)[0]
            read_start = division_coor[x][0]
            read_end = division_coor[x][1]

        chr_data = read_bam(SSCS_bam,
                            pair_dict=pair_dict,
                            read_dict=read_dict,
                            csn_pair_dict=csn_pair_dict,
                            tag_dict=tag_dict,
                            badRead_bam=badRead_bam,
                            duplex=True,
                            read_chr=read_chr,
                            read_start=read_start,
                            read_end=read_end
                            )

        read_dict = chr_data[0]
        tag_dict = chr_data[1]
        pair_dict = chr_data[2]
        csn_pair_dict = chr_data[3]

        counter += chr_data[4]
        unmapped += chr_data[5]
        unmapped_mate += chr_data[6]
        multiple_mapping += chr_data[7]

        # ===== Create consenus seq for reads =====
        for readPair in list(csn_pair_dict.keys()):
            if len(csn_pair_dict[readPair]) != 2:
                print('uh oh pairing problem!!! Only one read, mate missing')
                print(readPair)
                print(csn_pair_dict[readPair][0])
                print(read_dict[csn_pair_dict[readPair][0]][0])
                print(SSCS_bam.mate(read_dict[csn_pair_dict[readPair][0]][0]))
            else:
                for tag in csn_pair_dict[readPair]:
                    # === Determine tag of duplex pair read ===
                    ds = duplex_tag(tag)

                    # === Group duplex read pairs and create consensus ===
                    # Check presence of duplex pair
                    if ds not in duplex_dict.keys():
                        if tag in tag_dict and ds in tag_dict:
                            duplex_count += 1

                            # consensus seq
                            consensus_seq, qual_consensus = duplex_consensus(read_dict[tag][0], read_dict[ds][0])

                            # consensus duplex tag
                            dcs_query_name = dcs_consensus_tag(tag, ds)

                            dcs_read = create_aligned_segment([read_dict[tag][0], read_dict[ds][0]], consensus_seq,
                                                              qual_consensus, dcs_query_name)

                            duplex_dict[tag] += 1
                            # duplex_dict[tag] = dcs_read  # add duplex tag to dictionary to prevent making a duplex for the same sequences twice

                            DCS_bam.write(dcs_read)

                        else:
                            SSCS_singleton.write(read_dict[tag][0])
                            SSCS_singletons += 1

                        # Remove read from dictionary after writing
                        del read_dict[tag]

                # Remove key from dictionary after writing
                del csn_pair_dict[readPair]

    print(collections.Counter([i for i in tag_dict.values()]))  ## WHY ARE SOME VALUES 0 IN THE COUNTER?
    time_tracker.write('DCS: ')
    time_tracker.write(str((time.time() - start_time)/60) + '\n')

    summary_stats = '''
SSCS - Total reads: {} \n
SSCS - Unmapped reads: {} \n
SSCS - Secondary/Supplementary reads: {} \n
DCS reads: {} \n
SSCS singletons: {} \n'''.format(counter, unmapped, multiple_mappings, duplex_count, SSCS_singletons)
    stats.write(summary_stats)

    print(summary_stats)

    time_tracker.close()
    stats.close()
    DCS_bam.close()
    SSCS_singleton.close()

    return duplex_dict


###############################
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
