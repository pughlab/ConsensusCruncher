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

from consensus_helper_sc import *

###############################
##         Functions         ##
###############################

#def duplex_consensus(seq1, seq2, readlength):
    #'''(str, str, int) -> str
    #Return consensus of 2 sequences with N for variant bases.
    #'''    
    #consensus = ''
    #if seq1 == seq2:
        #consensus = seq1
    #else:
        #for i in range(readlength):
            #try:
                #if seq1[i] == seq2[i]:
                    #consensus += seq1[i]
                #else:
                    #consensus += 'N'
            #except:
                #consensus += 'N'
    
    #return consensus

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


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True, type = str)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output BAM file", required = True, type = str)
    args = parser.parse_args()
    
    start_time = time.time()    
   
    args.infile = str(args.infile)
    args.outfile = str(args.outfile)
 
    SSCS_bam = pysam.AlignmentFile(args.infile, "rb")    
    DCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = SSCS_bam)
    SSCS_singleton = pysam.AlignmentFile('{}.dcs.singleton.bam'.format(args.outfile.split('.dcs')[0]), "wb", template = SSCS_bam)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.dcs')[0]), "wb", template = SSCS_bam)

    stats = open('{}.stats.txt'.format(args.outfile.split('.dcs')[0]), 'a')
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.dcs')[0]), 'a')
    
    # setup fastq files
    # fastqFile1 = open('{}.dcs_R1.fastq.gz'.format(args.outfile.split('.dcs')[0]), 'w')
    # fastqFile2 = open('{}.dcs_R2.fastq.gz'.format(args.outfile.split('.dcs')[0]), 'w')
    
    # ===== Initialize dictionaries and convert SSCS data into dict=====    
    read_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    pair_dict = collections.defaultdict(list)
    csn_pair_dict = collections.defaultdict(list)
    
    unmapped = 0
    unmapped_flag = 0
    bad_reads = 0 # secondary/supplementary reads
    poor_mapq = 0
    counter = 0  
    SSCS_singletons = 0
    multiple_mappings = 0
    
    duplex_count = 0
    duplex_dict = collections.OrderedDict()    
    duplex_pair = collections.defaultdict(int)   
    
    # ===== Read data and create dictionaries =====
    chr_data = read_bam(SSCS_bam, 
                        pair_dict = pair_dict,
                        read_dict = read_dict,
                        tag_dict = tag_dict,
                        csn_pair_dict=csn_pair_dict,
                        badRead_bam = badRead_bam,
                        duplex = True)

    read_dict = chr_data[0]
    tag_dict = chr_data[1]
    pair_dict = chr_data[2]
    csn_pair_dict = chr_data[3]

    counter += chr_data[4]
    unmapped += chr_data[5]
    unmapped_flag += chr_data[6]
    bad_reads += chr_data[7]
    poor_mapq += chr_data[8]
    
    # ===== Create consenus seq for reads =====
    for readPair in list(csn_pair_dict.keys()):
        if len(csn_pair_dict[readPair]) != 2:
            print('uh oh pairing problem!!! Only one read, mate missing')
            print(readPair)
            print(csn_pair_dict[readPair])
        else:
            for tag in csn_pair_dict[readPair]:
                # === Determine tag of duplex pair read ===
                barcode = tag.split('_')[0]
                barcode_bases = int(len(barcode)/2)  # number of barcode bases, avoids complications if num bases change
                # duplex barcode is the reverse (e.g. AT|GC -> GC|AT [dup])
                duplex_barcode = barcode[barcode_bases:] + barcode[:barcode_bases]

                read_num = tag[-2:]
                if read_num == 'R1':
                    dup_read_num = 'R2'
                else:
                    dup_read_num = 'R1'

                ds = duplex_barcode + '_' + tag.split('_', 1)[1][:-2] + dup_read_num  # duplex tag

                # === Group duplex read pairs and create consensus ===
                # Check presence of duplex pair
                if ds not in duplex_dict.keys():
                    if tag_dict[tag] == 1 and tag_dict[ds] == 1:
                        print('ds')
                        print(tag)
                        print(ds)
                        duplex_count += 1
                        duplex_dict[ds] = tag

                        # consensus seq
                        consensus_seq, qual_consensus = duplex_consensus(read_dict[tag][0], read_tag[ds][0])

                        try:
                            dcs_read = create_aligned_segment([read], consensus_seq, qual_consensus)
                        except:
                            print(read.query_qualities)
                            print(duplex.query_qualities)
                            print(consensus_seq, qual_consensus)

                        dcs_read.query_name = "{}_{}_{}".format(min(barcode, pair_barcode), max(barcode, pair_barcode),
                                                                dcs_read.query_name.split('_', 1)[1].rsplit('_', 3)[0])

                        duplex_dict[tag] = dcs_read  # add duplex tag to dictionary to prevent making a duplex for the same sequences twice
                        duplex_pair[readPair] += 1

                        DCS_bam.write(dcs_read)

                    else:
                        SSCS_singleton.write(read_dict[tag][0])
                        SSCS_singletons += 1

    print(tag_dict)
    print(collections.Counter([i for i in tag_dict.values()]))  ## WHY ARE SOME VALUES 0 IN THE COUNTER?
    time_tracker.write('DCS: ')
    time_tracker.write(str((time.time() - start_time)/60) + '\n') 
    
    summary_stats = '''\n
=== DCS ===
Total reads: {} \n
Unmapped reads: {} \n
Unmapped flag reads: {} \n
Secondary/Supplementary reads: {} \n
Multiple Mappings: {} \n
DCS reads: {} \n
SSCS singletons: {} \n
    '''.format(counter, unmapped, unmapped_flag, bad_reads, multiple_mappings, duplex_count, SSCS_singletons)   
    stats.write(summary_stats)

    print(summary_stats)
    
    time_tracker.close()    
    stats.close()
    DCS_bam.close()
    SSCS_singleton.close()
    
    return duplex_dict


###############################
##           Main            ##
###############################
if __name__ == "__main__": 
    import time
    start_time = time.time()
    main()
