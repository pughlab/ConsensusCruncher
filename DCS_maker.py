#!/usr/bin/env python

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

from consensus_helper import *

###############################
##         Functions         ##
###############################

def duplex_consensus(seq1, seq2, readlength):
    '''(str, str, int) -> str
    Return consensus of 2 sequences with N for variant bases.
    '''    
    consensus = ''
    if seq1 == seq2:
        consensus = seq1
    else:
        for i in range(readlength):
            try:
                if seq1[i] == seq2[i]:
                    consensus += seq1[i]
                else:
                    consensus += 'N'
            except:
                consensus += 'N'
    
    return consensus


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output BAM file", required = True)
    args = parser.parse_args()
    
    start_time = time.time()    
    
    SSCS_bam = pysam.AlignmentFile(args.infile, "rb")    
    DCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = SSCS_bam)
    SSCS_singleton = pysam.AlignmentFile('{}.sscs.singleton.bam'.format(args.outfile.split('.dcs')[0]), "wb", template = SSCS_bam)
    stats = open('{}_stats.txt'.format(args.outfile.split('.dcs')[0]), 'a')
    time_tracker = open('{}_time_tracker.txt'.format(args.outfile.split('.dcs')[0]), 'a')
    
    
    bam_dicts = uid_dict(SSCS_bam)
    SSCS_dict = bam_dicts[0]
    
    duplex_count = 0
    duplex_dict = collections.OrderedDict()
    
    for i in SSCS_dict.keys():
        #print(i)
        barcode = i.split('_')[0]
        barcode_bases = int(len(barcode)/2)
        pair_barcode = barcode[barcode_bases:] + barcode[:barcode_bases]
        
        read_num = i[-2:]
        if read_num == 'R1':
            read_num = 'R2'
        else:
            read_num = 'R1'
                
        ds = pair_barcode + '_' + i.split('_', 1)[1][:-2] + read_num
        
        if ds in SSCS_dict.keys() and ds not in duplex_dict.keys():
            duplex_count += 1
            duplex_dict[ds] = i
            
            read = SSCS_dict[i][0]
            duplex = SSCS_dict[ds][0]
            #print(i)
            #print(read)
            #print(ds)
            #print(duplex)
            dsc = duplex_consensus(read.query_alignment_sequence, duplex.query_alignment_sequence, read.query_alignment_length)  
            
            dsc_read = create_aligned_segment([read, duplex], dsc, read.query_alignment_qualities)
            dsc_read.query_name = "{}|{}|{}\t".format(dsc_read.query_name.split('|')[0], min(barcode, pair_barcode), max(barcode, pair_barcode)) # Add both barcodes to header in alphabetical order
            
            duplex_dict[i] = dsc_read # when i is the duplexed pair, it'll find that i is already in dictionary
            
            DCS_bam.write(dsc_read)   
        else:
            SSCS_singleton.write(SSCS_dict[i][0])


    time_tracker.write('DCS: ')
    time_tracker.write(str((time.time() - start_time)/60) + '\n') 
    
    summary_stats = '''DCS reads: {}'''.format(duplex_count)   
    stats.write(summary_stats)
    
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
    print((time.time() - start_time)/60)
    

