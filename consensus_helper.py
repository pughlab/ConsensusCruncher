#!/usr/bin/env python

###############################################################
#
#                          Consensus Helper
#
# Author: Nina Wang
# Last Modified: May 18, 2016
# Date Created: Mar 24, 2016
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

###############################
##         Functions         ##
###############################

def uid_dict(bamfile, read_chr = None, read_start = None, read_end = None):
    '''(bamfile object) -> dict, dict, int, int
    
    Input: bamfile object created from pysam.AlignmentFile
    
    Output:
    1) bam_dict: dictionary of bamfile reads
                 - Key: barcode_chr_startR1_startR2_strand_ReadNum
                 - Value: bamfile read
    2) tag_dict: integer dictionary of number of reads in each read family
    3) counter: total number of reads
    4) bad_reads: number of reads that not properly mapped (lack cigar string)
    '''
    bam_dict = collections.OrderedDict() #dict that remembers order entries were added
    tag_dict = collections.defaultdict(int) #dict tracking int
    
    if read_chr == None:
        bamLines = bamfile.fetch(until_eof = True)
    else:
        bamLines = bamfile.fetch(read_chr, read_start, read_end)
    
    bad_reads = 0
    counter = 0
    
    for line in bamLines:
        #print(line)
        counter += 1
        
        strand = 'fwd'
        if line.is_reverse:
            strand = 'rev'
        
        read = 'R1'
        if line.is_read2:
            read = 'R2'
            
        try:
            tag = '{}_{}_{}_{}_{}_{}'.format(line.qname.split("|")[1], # mol barcode
                                          line.reference_name, # chr num
                                          line.reference_start, # Start R1 (0-based)
                                          line.next_reference_start, # Start R2
                                          strand, # strand direction
                                          read # Read num
                                          )        
            
            # Raise error if cigarstring is empty indicating bad read
            if 'I' in line.cigarstring:
                line.cigarstring

            tag_dict[tag] += 1     
            
            if tag not in bam_dict:
                bam_dict[tag] =[line]
        
            else:
                bam_dict[tag].append(line)             
            
        except:
            # Bad reads won't have cigar or MD 
            bad_reads += 1

    return bam_dict, tag_dict, counter, bad_reads


def read_mode(field, bam_reads):
    '''(str, lst) -> str
    Return mode (most common occurrence) of specified field (e.g. cigarstring, flag, mapping quality, template_length).
    '''
    field = 'i.{}'.format(field)
    # Rank by number of occurrences 
    field_lst = collections.Counter(eval(field) for i in bam_reads).most_common() 
    # Take max occurrences
    common_field_lst = [i for i, j in field_lst if j == field_lst[0][1]] 
    # Randomly select max if there's multiple   
    common_field = common_field_lst[randint(0, len(common_field_lst)-1)] 
    
    return common_field


def create_aligned_segment(bam_reads, sscs, sscs_qual):
    '''(list, str, list) -> pysam object
    Return pysam object with new consensus seq given list of bam reads.
    
    '''
    # Find most common cigar seq
    common_cigar = read_mode('cigarstring', bam_reads)
    
    # Find first read with most common cigar and set as template read for SSCS
    template_index = [i.cigarstring for i in bam_reads].index(common_cigar)    
    template_read = bam_reads[template_index]
    #print(template_read)
    
    # Create bam read based on template read
    SSCS_read = pysam.AlignedSegment()
    SSCS_read.query_name = template_read.query_name
    SSCS_read.flag = read_mode('flag', bam_reads) # Take most common flag
    SSCS_read.query_sequence = sscs
    SSCS_read.reference_id = template_read.reference_id
    SSCS_read.reference_start = template_read.reference_start
    SSCS_read.mapping_quality = read_mode('mapping_quality', bam_reads) # Most common mapping quality
    SSCS_read.cigar = template_read.cigar
    SSCS_read.next_reference_id= template_read.next_reference_id
    SSCS_read.next_reference_start = template_read.next_reference_start  
    SSCS_read.template_length = read_mode('template_length', bam_reads)
    SSCS_read.query_qualities = sscs_qual
    
    SSCS_read.set_tag('RG', read_mode("get_tag('RG')", bam_reads))
    #SSCS_read.set_tag('MD', )
    
    
    # --NOTE: Tags currently disabled as it gives errors when bams are loaded into IGV 
    #print(read_mode("get_tag('RG')", bam_reads))
    #SSCS_read.tags = read_mode("get_tag('RG')", bam_reads)
    #tag_index = [x for x, y in enumerate(SSCS_read.tags) if y[0] == 'MD'][0]
    #SSCS_read.tags[tag_index] = read_mode("get_tag('MD')", bam_reads)
    #SSCS_read.tags[tag_index] = read_mode("get_tag('RG')", bam_reads)
    
    return SSCS_read




#filename = 'MEM-001_KRAS.sscs.bam'
#SSCS_bam = pysam.AlignmentFile(filename, "rb")

