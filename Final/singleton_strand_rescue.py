#!/usr/bin/env python3

###############################################################
#
#                      Singleton Rescue
#
#  Author: Nina Wang
#  Date Created: July 5, 2016
###############################################################
#  Function:
#  Written for Python 3.5.1
#
#  uid_dict(bamfile): create bam_read dictionaries and counts number of reads
#  read_mode(field, bam_reads): return most common occurrence in a specified field of read (e.g. cigar, etc)
#  create_aligned_segment(bam_reads, sscs, sscs_qual): Create consensus bam read
#
#  Inputs:
#  1. A position-sorted paired-end BAM file containing reads with a duplex tag
#    in the header.
#
#  Description:
#  Helper functions for single strand and duplex consensus making
#
###############################################################

import pysam # Need to install
import collections
from argparse import ArgumentParser
import time
import math
import os
import inspect

from consensus_helper import *
#  from SSCS_maker import consensus_maker


###############################
##         Functions         ##
###############################
def duplex_consensus(read1, read2):
    '''(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> pysam.calignedsegment.AlignedSegment

    Return consensus of 2 reads with N for variant bases.
    '''
    consensus_seq = ''
    qual_consensus = []

    ### NOTE: if one of the reads is SSCS, maybe we should take the base of that instead? AKA just asign the singleton read the seq of the SSCS
    ### My concern is when the SSCS is formed from 2/3 reads...

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


def duplex_tag(tag):
    '''(str) -> str
    Return duplex tag of given tag.

    Example:
    tag = TTCT_24_1584331_24_1584409_98M_98M_fwd_R1
    duplex = CTTT_24_1584331_24_1584409_98M_98M_fwd_R2
    '''
    barcode = tag.split('_')[0]
    barcode_bases = int(len(barcode) / 2)
    pair_barcode = barcode[barcode_bases:] + barcode[:barcode_bases]  # get barcode of mate pair read (opposite strand)

    read_num = tag[-2:]
    if read_num == 'R1':
        read_num = 'R2'
    else:
        read_num = 'R1'

    strand_ori = tag.split('_')[7]
    if strand_ori == 'pos':
        duplex_ori = 'neg'
    else:
        duplex_ori = 'pos'

    # pair read is other read on opposite strand (e.g. read = ACGT_fwd_R1, duplex = GTAC_fwd_R2)
    duplex = pair_barcode + '_' + tag.split('_', 1)[1][:-2] + read_num
    duplex = duplex.replace(strand_ori, duplex_ori)

    return duplex


def strand_rescue(read_tag, duplex_tag, query_name, singleton_dict, sscs_dict=None):
    """(str, str, dict, dict) -> Pysam.AlignedSegment

    Return 'rescued' singleton read using compliment read from opposite strand (either found in SSCS or singleton).

    Quality score calculated from singleton and complimentary read. Read template retained from singleton being rescued.
    """
    read = singleton_dict[read_tag][0]

    # If SSCS bamfile provided, rescue with SSCS
    # try:
    if sscs_dict == None:
        compliment_read = singleton_dict[duplex_tag][0]
    else:
        compliment_read = sscs_dict[duplex_tag][0]
    # except KeyError:
        # If duplex not found in dictionary, output None
        # return None

    dcs = duplex_consensus(read, compliment_read)

    dcs_read = create_aligned_segment([read], dcs[0], dcs[1])
    dcs_read.query_name = query_name

    return dcs_read


def main():
    '''Singleton rescue:
    - First rescue with SSCS bam
    - Rescue remaining singletons with singleton bam
    '''
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--singleton", action = "store", dest="singleton", help="input singleton BAM file",
                        required = True, type = str)
    # parser.add_argument("--sscs", action = "store", dest="sscs", help="input singleton BAM file", required = False)
    #parser.add_argument("--rescue_outfile", action = "store", dest="rescue_outfile", help="output BAM file", required = True)
    parser.add_argument("--bedfile", action = "store", dest="bedfile", help="input bedfile", required = False)
    args = parser.parse_args()

    start_time = time.time()

    ######################
    ##      SETUP       ##
    ######################
    # Read input bams as pysam objects
    singleton_bam = pysam.AlignmentFile(args.singleton, "rb")
    sscs_bam = pysam.AlignmentFile('{}.sscs.sort.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), "rb")

    # Setup output bams
    sscs_rescue_bam = pysam.AlignmentFile('{}.sscs.rescue.bam'.format(args.singleton.split('.singleton.sort.bam')[0]),
                                          'wb',
                                          template = singleton_bam)
    singleton_rescue_bam = pysam.AlignmentFile('{}.singleton.rescue.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), 'wb',
                                          template=singleton_bam)
    remaining_rescue_bam = pysam.AlignmentFile(
        '{}.rescue.remaining.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), 'wb',
        template=singleton_bam)

    badRead_bam = pysam.AlignmentFile('{}.singleton.badReads.bam'.format(args.singleton.split('.singleton.sort.bam')[0]),
                                      "wb", template = singleton_bam)
    stats = open('{}.rescue_stats.txt'.format(args.singleton.split('.singleton.sort.bam')[0]), 'w')

    # set up time tracker
    time_tracker = open('{}.time.tracker.txt'.format(args.singleton.split('.singleton.sort.bam')[0]), 'a')

    # Initialize dictionaries
    singleton_dict = collections.OrderedDict() # dict that remembers order of entries
    singleton_tag = collections.defaultdict(int)
    singleton_pair = collections.OrderedDict()
    singleton_csn_pair = collections.defaultdict(list)

    sscs_dict = collections.OrderedDict() # dict that remembers order of entries
    sscs_tag = collections.defaultdict(int)
    sscs_pair = collections.OrderedDict()
    sscs_csn_pair = collections.defaultdict(list)

    rescue_dict = collections.OrderedDict()

    # Initialize counters
    singleton_counter = 0 # this differs from 'counter' (which represents number of reads in dict as we're using
    # pair_dict, meaning only pairs are retained) => consensus_helper script may not be sufficient at making read pairs?
    singleton_unmapped = 0
    singleton_bad_reads = 0

    sscs_counter = 0
    sscs_unmapped = 0
    sscs_bad_reads = 0

    sscs_dup_rescue = 0
    singleton_dup_rescue = 0
    singleton_remaining = 0

    counter = 0

    #######################
    ##  SPLIT BY REGION  ##
    #######################

    chrm = [x['SN'] for x in singleton_bam.header['SQ']]
    chr_len = [x['LN'] for x in singleton_bam.header['SQ']]

    if 'args.bedfile' in locals():
        division_coor = coor_separator(chrm, chr_len, args.bedfile)
    else:
        division_coor = coor_separator(chrm, chr_len)

    for x in division_coor.keys():
        # Create dictionaries from bamfiles
        singleton = read_bam(singleton_bam,
                             pair_dict = singleton_pair,
                             read_dict = singleton_dict, # keeps track of paired tags
                             tag_dict = singleton_tag,
                             csn_pair_dict=csn_pair_dict,
                             badRead_bam=badRead_bam,
                             read_chr=x.rsplit('_', 1)[0],
                             read_start=division_coor[x][0],
                             read_end=division_coor[x][1]
                             )

        singleton_dict = singleton[0]
        singleton_tag = singleton[1]
        singleton_pair = singleton[2]
        singleton_csn_pair = singleton[3]

        singleton_counter += singleton[4]
        singleton_unmapped += singleton[5]
        singleton_bad_reads += singleton[6]

        sscs = read_bam(sscs_bam,
                        pair_dict = sscs_pair,
                        read_dict = sscs_dict,  # keeps track of paired tags
                        tag_dict = sscs_tag,
                        csn_pair_dict=csn_pair_dict,
                        badRead_bam=badRead_bam,
                        duplex = True,
                        read_chr=x.rsplit('_', 1)[0],
                        read_start=division_coor[x][0],
                        read_end=division_coor[x][1]
                        )

        sscs_dict = sscs[0]
        sscs_tag = sscs[1]
        sscs_pair = sscs[2]
        sscs_csn_pair = sscs[3]

        sscs_counter += sscs[4]
        sscs_unmapped += sscs[5]
        sscs_bad_reads += sscs[6]


        ######################
        ##      RESCUE      ##
        ######################
        for readPair in list(singleton_pair.keys()):
            # print(readPair)
            # print(singleton_pair[readPair])
            # return(readPair)

            for tag in singleton_pair[readPair]:
                counter += 1
                # print(readPair)
                # print(tag)
                # Check to see if singleton can be rescued by SSCS, then by singletons. If not, add to 'remaining' rescue bamfile
                duplex = duplex_tag(tag)
                query_name = readPair + ':' + str(singleton_tag[tag])

                # print(duplex)
                # return 'hi'

                # if 'CTCT_11_108143378_11_108143521_147M_147M_fwd_R2' in singleton_dict:
                #     print(duplex in sscs_dict)
                #     print(duplex in singleton_dict)

                # 1) SSCS duplex strand rescue -> check for duplex sequence matching singleton in SSCS bam
                if duplex in sscs_dict.keys():
                    rescue_read = strand_rescue(tag, duplex, query_name, singleton_dict, sscs_dict=sscs_dict)
                    sscs_dup_rescue += 1
                    sscs_rescue_bam.write(rescue_read)

                    # rescue_dict[tag] = duplex

                    del sscs_dict[duplex]
                    del singleton_dict[tag]

                # 2) Singleton duplex strand rescue -> check for duplex sequence matching singleton in singletons bam
                elif duplex in singleton_dict.keys():
                    rescue_read = strand_rescue(tag, duplex, query_name, singleton_dict)
                    singleton_dup_rescue += 1
                    singleton_rescue_bam.write(rescue_read)
                    rescue_dict[tag] = duplex

                    if duplex in rescue_dict.keys():
                        # print(tag)
                        # print(duplex)
                        # print(rescue_dict)
                        # print(tag in singleton_dict)
                        # print(duplex in singleton_dict)
                        # return 'hi'
                        # delete tags from dict if both singletons are rescued in singleton-singleton
                        # (aka duplex strand is already in dict)
                        del singleton_dict[tag]
                        del singleton_dict[duplex]

                # 3) Singleton written to remaining bam if neither SSCS or Singleton duplex rescue was possible
                else:
                    remaining_rescue_bam.write(singleton_dict[tag][0])
                    singleton_remaining += 1
                    del singleton_dict[tag]

            del singleton_pair[readPair]

        # Update time tracker
        time_diff = str((time.time() - start_time)/60)
        print(time_diff)
        try:
            time_tracker.write(x + ': ')
            # time_tracker.write(time_diff + '\n')
        except:
            continue

    ######################
    ##      SUMMARY     ##
    ######################
    sscs_rescue_frac = (sscs_dup_rescue/singleton_counter) * 100
    singleton_rescue_frac = (singleton_dup_rescue/singleton_counter) * 100

    summary_stats = '''Total singletons: {} \n
SSCS strand rescued singletons: {} \n
% SSCS rescue: {} \n
Singleton strand rescued singletons: {} \n
% singleton rescue: {} \n
Singletons remaining (not rescued): {} \n
'''.format(counter, sscs_dup_rescue, sscs_rescue_frac, singleton_dup_rescue, singleton_rescue_frac, singleton_remaining)

    print(singleton_counter)
    print(sscs_dup_rescue)
    print(singleton_dup_rescue)
    print(singleton_remaining)
    print(counter)

    stats.write(summary_stats)

    # Close files
    singleton_bam.close()
    sscs_bam.close()
    sscs_rescue_bam.close()
    singleton_rescue_bam.close()
    remaining_rescue_bam.close()
    stats.close()


###############################
##           Main            ##
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)