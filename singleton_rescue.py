#!/usr/bin/env python

###############################################################
#
#                      Singleton Rescue
#
# Author: Nina Wang
# Date Created: July 5, 2016
###############################################################
# Function: 
# Written for Python 3.5.1
# 
# uid_dict(bamfile): create bam_read dictionaries and counts number of reads
# read_mode(field, bam_reads): return most common occurrence in a specified field of read (e.g. cigar, etc)
# create_aligned_segment(bam_reads, sscs, sscs_qual): Create consensus bam read
#
# Inputs: 
# 1. A position-sorted paired-end BAM file containing reads with a duplex tag 
#    in the header.  
#
# Description:
# Helper functions for single strand and duplex consensus making
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


def rescue_duplex(read_tag, duplex_tag, singleton_dict, bamfile = None):
    '''(str, str, dict, dict) -> Pysam.AlignedSegment
    
    Return consensus read from singleton + duplex read.
    - Duplex read: Matching sequence on the other strand
    
    Quality scores: SSCS quality scores or average of two singletons
    
    If not additional bamfile provided, singleton_dict will bee used for read rescue.
    '''
    
    read = singleton_dict[read_tag][0]
    
    # If bamfile provided, perform SSCS singleton rescue 
    if bamfile == None:
        duplex_read = singleton_dict[duplex_tag][0]
    else:
        duplex_read = bamfile[duplex_tag][0]
        
    dcs = duplex_consensus(read.query_alignment_sequence, duplex_read.query_alignment_sequence, read.query_alignment_length)

    # If consensus has >30% N's, toss read
    if dcs.count('N')/len(dcs) > 0.3:
        return None
    
    # If SSCS used to rescue read, create segment based on 
    if bamfile == None:
        dsc_read = create_aligned_segment([read], dcs, read.query_alignment_qualities)
    else:
        dsc_read = create_aligned_segment([read, duplex_read], dcs, read.query_alignment_qualities)
        

    return dsc_read


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--infile", action = "store", dest="infile", help="input singleton BAM file", required = True)
    #parser.add_argument("--rescue_outfile", action = "store", dest="rescue_outfile", help="output BAM file", required = True)
    args = parser.parse_args()
    
    ### Later include flags to turn on/off certain rescues ###
    
    
    start_time = time.time()
    
    singleton_bam = pysam.AlignmentFile(args.infile, "rb") 
    sscs_bam = pysam.AlignmentFile('{}.sscs.bam'.format(args.infile.split('.singleton.sort.bam')[0]), "rb")
    rescue_bam = pysam.AlignmentFile('{}.rescue.bam'.format(args.infile.split('.singleton.sort.bam')[0]), "wb", template = sscs_bam)
    stats = open('{}_stats.txt'.format(args.infile.split('.singleton.sort.bam')[0]), 'a')
    time_tracker = open('{}_time_tracker.txt'.format(args.infile.split('.singleton.sort.bam')[0]), 'a')
    
    singleton_dict = uid_dict(singleton_bam)[0]
    sscs_dict = uid_dict(sscs_bam)[0]
    
    
    sscs_dup_rescue = 0
    singleton_dup_rescue = 0
    
    rescue_dict = collections.OrderedDict()
    
    for i in singleton_dict.keys():
        barcode = i.split('_')[0]
        barcode_bases = int(len(barcode)/2)
        pair_barcode = barcode[barcode_bases:] + barcode[:barcode_bases] # get barcode of mate pair read (opposite strand)
        
        read_num = i[-2:]
        if read_num == 'R1':
            read_num = 'R2'
        else:
            read_num = 'R1'
        
        duplex = pair_barcode + '_' + i.split('_', 1)[1][:-2] + read_num # pair read is other read on opposite strand
        
        
        # Check if there's duplicate singletons or singletons already rescued
        if read not in rescue_dict.keys():
        
            # 1) Duplex rescue: Singleton + SSCS -> Check for singleton matching duplex sequence in SSCS bam
            if duplex in sscs_dict.keys():            
                dcs_read = rescue_duplex(i, duplex, singleton_dict, sscs_dict)           
                # Check if read passes Ncutoff (0.3)
                if dcs_read != None:
                    sscs_dup_rescue += 1
                    #dcs_read.query_name = "{}|{}|{}\t".format(dcs_read.query_name.split('|')[0], min(barcode, pair_barcode), max(barcode, pair_barcode)) # Add both barcodes to header in alphabetical order [NOT CURRENTLY IMPLEMENTED - as rescued reads still go through duplex consensus making]
                    rescue_bam.write(dcs_read)
                    
            # 2) Duplex rescue: Singleton + Singleton -> Check for singleton matching duplex sequence in Singleton bam -> WOULD YOU RESCUE BOTH SINGLETONS THEN??
            elif duplex in singleton_dict.keys() and duplex not in rescue_dict.keys():            
                dcs_read = rescue_duplex(i, duplex, singleton_dict)
                
                if dcs_read != None:
                    singleton_dup_rescue += 1                
                    #dcs_read.query_name = "{}|{}|{}\t".format(dcs_read.query_name.split('|')[0], min(barcode, pair_barcode), max(barcode, pair_barcode)) # Add both barcodes to header in alphabetical order                 
                    rescue_bam.write(dcs_read)
            
            # 3) Non-duplex rescue: Any SSCS mapping to region 
            
            #else:

# rescued bams should be based on the singleton sequence, but with the duplex consensus 
# duplex consensus should have filters (Ncutoff) to reflect 


    time_tracker.write('Singleton Rescue: ')
    time_tracker.write(str((time.time() - start_time)/60) + '\n')

    
    summary_stats='''SSCS duplex rescue: {} \n
Singleton duplex rescue: {} \n'''.format(sscs_dup_rescue, singleton_dup_rescue)
    
    stats.write(summary_stats)
    
    
    singleton_bam.close()
    sscs_bam.close()
    rescue_bam.close()
    stats.close()
    time_tracker.close()
    

###############################
##           Main            ##
###############################
if __name__ == "__main__": 
    import time
    start_time = time.time()
    main()  
    print((time.time() - start_time)/60)  
