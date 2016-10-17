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

from consensus_helper_sc import *
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

    duplex = pair_barcode + '_' + tag.split('_', 1)[1][
                                  :-2] + read_num  # pair read is other read on opposite strand (e.g. read = ACGT_fwd_R1, duplex = GTAC_fwd_R2)
    return duplex


def strand_rescue(read_tag, duplex_tag, singleton_dict, sscs_dict=None):
    """(str, str, dict, dict) -> Pysam.AlignedSegment

    Return 'rescued' singleton read using compliment read from opposite strand (either found in SSCS or singleton).

    Quality score calculated from singleton and complimentary read. Read template retained from singleton being rescued.
    """
    read = singleton_dict[read_tag][0]

    # If SSCS bamfile provided, rescue with SSCS
    try:
        if sscs_dict == None:
            compliment_read = singleton_dict[duplex_tag][0]
        else:
            compliment_read = sscs_dict[duplex_tag][0]
    except KeyError:
        # If duplex not found in dictionary, output None
        return None

    dcs = duplex_consensus(read, compliment_read)

    dcs_read = create_aligned_segment([read], dcs[0], dcs[1])

    return dcs_read


def cytoband_pos(chr_lst, chr_len, bedfile = None):
    '''(list, list) -> list
    Return list of int indicating chromosomal arm positions given a list of chromosomes and their lengths.

    ChrM not divided by arms.

    Chrm arm positions are used to separate bam file reads into more manageable chunks, so dictionaries don't take up too much memory.

    Input:
    - chr_lst
    ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    - chr_len
    [16571, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]

    '''
    cytoband_coor = collections.OrderedDict()

    if 'chrM' in chr_lst:
        chr_arm_coor['chrM'] = (0, chr_len[chr_lst.index('chrM')])

    if bedfile == None:
        filepath = os.path.abspath(inspect.getfile(inspect.currentframe())).rsplit('/', 1)[0]
        bedfile = filepath + '/cytoBand.txt'

    ### Should write script to incorporate any chr not found in cytoband file....!!!

    with open(bedfile) as f:
        next(f) # Skip header
        for line in f:
            chr_arm = line.split('\t')
            chr_key = '{}_{}'.format(chr_arm[0], chr_arm[3])
            start = int(chr_arm[1]) # python is 0-based (start is usually 1)
            end = int(chr_arm[2])
            chr_val = (start, end)

            cytoband_coor[chr_key] = chr_val

    #if 'chrHPV16_gi_333031' in chr_lst:
        #chr_arm_coor['chrHPV16_gi_333031'] = (0, chr_len[chr_lst.index('chrHPV16_gi_333031')])

    return cytoband_coor


def main():
    '''Singleton rescue:
    - First rescue with SSCS bam
    - Rescue remaining singletons with singleton bam
    '''
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--singleton", action = "store", dest="singleton", help="input singleton BAM file", required = True)
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
    sscs_bam = pysam.AlignmentFile('{}.gthan2.sscs.sort.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), "rb")

    # Setup output bams
    sscs_rescue_bam = pysam.AlignmentFile('{}.sscs.rescue.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), 'wb',
                                     template = singleton_bam)
    singleton_rescue_bam = pysam.AlignmentFile('{}.singleton.rescue.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), 'wb',
                                          template=singleton_bam)
    remaining_rescue_bam = pysam.AlignmentFile(
        '{}.rescue.remaining.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), 'wb',
        template=singleton_bam)
    badRead_bam = pysam.AlignmentFile('{}.singleton.badReads.bam'.format(args.singleton.split('.singleton.sort.bam')[0]), "wb", template = singleton_bam)
    stats = open('{}.rescue_stats.txt'.format(args.singleton.split('.singleton.sort.bam')[0]), 'a')

    # Initialize dictionaries
    read_dict = collections.OrderedDict() # dict that remembers order of entries
    pair_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    read_pair_dict = collections.defaultdict(list)

    rescue_dict = collections.OrderedDict()

    # Initialize counters
    singleton_counter = 0
    singleton_unmapped = 0
    singleton_unmapped_flag = 0
    singleton_bad_reads = 0

    sscs_counter = 0
    sscs_unmapped = 0
    sscs_unmapped_flag = 0
    sscs_bad_reads = 0

    sscs_dup_rescue = 0
    singleton_dup_rescue = 0
    singleton_remaining = 0

    #######################
    ##  SPLIT BY REGION  ##
    #######################

    chrm = [x['SN'] for x in bamfile.header['SQ']]
    chr_len = [x['LN'] for x in bamfile.header['SQ']]

    if 'args.bedfile' in locals():
        cytoband_coor = cytoband_pos(chrm, chr_len, args.bedfile)
    else:
        cytoband_coor = cytoband_pos(chrm, chr_len)

    for x in cytoband_coor.keys():
        # Create dictionaries from bamfiles
        singleton = read_bam(singleton_bam,
                             pair_dict = pair_dict,
                             read_dict = read_dict, # keeps track of paired tags
                             tag_dict = tag_dict,
                             read_pair_dict = read_pair_dict, # keep track of paired reads (using query name),
                             #  reads removed from dict as they are added to pair_dict
                             badRead_bam= badRead_bam,
                             read_chr=x.rsplit('_', 1)[0],
                             read_start=chr_arm_coor[x][0],
                             read_end=chr_arm_coor[x][1]
                             )

        singleton_dict = singleton[0]
        singleton_tag = singleton[1]
        singleton_pair = singleton[2]
        singleton_read_pair = singleton[3]

        singleton_counter += singleton[4]
        singleton_unmapped += singleton[5]
        singleton_unmapped_flag += singleton[6]
        singleton_bad_reads += singleton[7]


        read_dict = collections.OrderedDict()  # dict that remembers order of entries
        pair_dict = collections.OrderedDict()
        tag_dict = collections.defaultdict(int)
        read_pair_dict = collections.defaultdict(list)

        sscs = read_bam(sscs_bam,
                        pair_dict = pair_dict,
                        read_dict = read_dict,  # keeps track of paired tags
                        tag_dict = tag_dict,
                        read_pair_dict = read_pair_dict,  # keep track of paired reads (using query name),
                        # reads removed from dict as they are added to pair_dict
                        badRead_bam = badRead_bam,
                        duplex = True,
                        read_chr=x.rsplit('_', 1)[0],
                        read_start=chr_arm_coor[x][0],
                        read_end=chr_arm_coor[x][1]
                        )

        sscs_dict = sscs[0]
        sscs_tag = sscs[1]
        sscs_pair = sscs[2]
        sscs_read_pair = sscs[3]

        sscs_counter += sscs[4]
        sscs_unmapped += sscs[5]
        sscs_unmapped_flag += sscs[6]
        sscs_bad_reads += sscs[7]


        ######################
        ##      RESCUE      ##
        ######################

        for tag in singleton_dict.keys():
            # Check to see if singleton can be rescued by SSCS, then by singletons. If not, add to 'remaining' rescue bamfile
            duplex = duplex_tag(tag)
            # print(tag)
            # print(duplex)
            #
            # print(duplex in sscs_dict.keys())
            # print(duplex in singleton_dict.keys())
            #
            # # print(sscs_dict.keys())
            # # print(singleton_dict.keys())
            # print(sscs_dict.keys() == singleton_dict.keys())
            # break

            # 1) SSCS duplex strand rescue -> check for duplex sequence matching singleton in SSCS bam
            if duplex in sscs_dict.keys():
                rescue_read = strand_rescue(tag, duplex, singleton_dict, sscs_dict = sscs_dict)
                if rescue_read != None:
                    sscs_dup_rescue += 1
                    sscs_rescue_bam.write(rescue_read)
                else:
                    singleton_remaining += 1
                    remaining_rescue_bam.write(singleton_dict[tag][0])
            # 2) Singleton duplex strand rescue -> check for duplex sequence matching singleton in singletons bam
            elif duplex in singleton_dict.keys():
                rescue_read = strand_rescue(tag, duplex, singleton_dict)
                if rescue_read != None:
                    singleton_dup_rescue += 1
                    singleton_rescue_bam.write(rescue_read)
                else:
                    singleton_remaining += 1
                    remaining_rescue_bam.write(singleton_dict[tag][0])
            # 3) Singleton written to remaining bam if neither SSCS or Singleton duplex rescue was possible
            else:
                singleton_remaining += 1
                remaining_rescue_bam.write(singleton_dict[tag][0])


    ######################
    ##      SUMMARY     ##
    ######################
    sscs_rescue_frac = sscs_dup_rescue/singleton_counter
    singleton_rescue_frac = singleton_dup_rescue/singleton_counter

    summary_stats = '''Total singletons: {} \n
SSCS strand rescued singletons: {} \n
% SSCS rescue: {} \n
Singleton strand rescued singletons: {} \n
% singleton rescue: {} \n
Singletons remaining (not rescued): {} \n
'''.format(singleton_counter, sscs_dup_rescue, sscs_rescue_frac, singleton_dup_rescue, singleton_rescue_frac, singleton_remaining)

    print(singleton_counter)
    print(sscs_dup_rescue)
    print(singleton_dup_rescue)

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