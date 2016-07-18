#!/usr/bin/env python

###############################################################
#
#                        Median SSCS Maker
#
# Author: Nina Wang
# Date Created: Jun 23, 2016
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
# 2. A singleton BAM file containing all read families with single reads.
# 3. A tag family size distribution plot (x-axis: family size, y-axis: number of reads).
# 4. A text file containing summary statistics (Total reads, SSCS reads, 
#    singletons, rescued reads)
#    - Read family: reads that share the same molecular barcode, chr, and start
#                   coordinates for Read1 and Read2
#    - Singleton: a read family containing only one member (a single read)
#
# usage: SSCS_maker.py [--cutoff CUTOFF] [--Ncutoff NCUTOFF] [--infile INFILE] 
#                      [--outfile OUTFILE]
# optional arguments:
# --cutoff CUTOFF     Percentage of nucleotides at a given position in a 
#                     sequence required to be identical for a consensus [0.7]
#                        Example (--cutoff = 0.7):
#                           Four reads (readlength = 10) are as follows:
#                              Read 1: ACTGATACTT
#                              Read 2: ACTGAAACCT
#                              Read 3: ACTGATACCT
#                              Read 4: ACTGATACTT
#                           The resulting SSCS is: ACTGATACNT    
#
# --Ncutoff NCUTOFF   Percentage of Ns allowed in a consensus sequence [0.3]
#                        Example (--ncutoff = 0.3):
#                           SSCS 1: ACGTGANCTAGTNCTNTACC
#                           SSCS 2: GATCTAGTNCATGACCGATA
#                        SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, 
#                        while SSCS 1 does not with 3/20 = 15% Ns.
#
# --infile INFILE     input BAM file
# --outfile OUTFILE   output BAM file
#
# python3 SSCS_maker.py --Ncutoff 0.3 --cutoff 0.7 --infile MEM-001_KRAS.bam --outfile MEM-001_KRAS.sscs.bam
#
###############################################################

import pysam # Need to install
import collections
import re
import array
from random import randint
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import statistics
import math

import inspect
import os

from consensus_helper import *

###############################
##         Functions         ##
###############################

def mismatch_pos(cigar, mismatch_tag):
    '''(list, str) -> lst
    Return 0-based index of mismatch positions in sequence (including insertions,
    ignoring deletions).
    
    E.g. mismatch_pos([(0, 19), (2, 1), (0, 79)], '19^A8G70')
    [27]
    First mismatch "G' at position 19 + 8
    
    Test cases:
    >>> mismatch_pos([(0, 98)], '31A0T22T41G') 
    [31, 32, 55, 97]
    >>> mismatch_pos([(0, 98)], '31A0T22T42')
    [31, 32, 55]
    >>> mismatch_pos([(0, 98)], 'G30A0T65')
    [0, 31, 32]
    
    In/del examples (1 = insertion, 2 = deletion in cigar):
    >>> mismatch_pos('19M1D26M1D28M25S', '19^A8G5T11^T8C0T12T5') # two deletions
    [27, 33, 53, 54, 67]
    >>> mismatch_pos([(0, 19), (2, 1), (0, 54), (4, 25)], '19^A8G5T10T8C0T12T5')
    [27, 33, 44, 53, 54, 67]
    >>> mismatch_pos([(0, 28), (1, 1), (0, 69)], '19T74G2') # one insertion 
    [19, 28, 94]
    >>> mismatch_pos([(0, 26), (1, 1), (0, 5), (2, 1), (0, 66)], '21C4G0G3^C14C6G5T33T4')
    [21, 26, 27, 45, 52, 58, 92]
    >>> mismatch_pos([(0, 74), (2, 2), (0, 3), (1, 2), (0, 19)], '70T1A1^GC22') # cases when insertions and deletions are in the same read
    [70, 72, 77]    
    >>> mismatch_pos([(0, 74), (2, 2), (0, 3), (2, 1), (1, 2), (0, 19)], '70T1A1^GC22')
    [70, 72, 77]
    >>> mismatch_pos('3M1I48M15D8M1D38M', '8G4C7G1G23C3^GAATTAAGAGAAGCA8^G38')
    >>> mismatch_pos([(0, 3), (1, 1), (0, 48), (2, 15), (0, 8), (2, 1), (0, 38)], '8G4C7G1G23C3^GAATTAAGAGAAGCA8^G38')
    [3, 8, 13, 21, 23, 47]
    
    Hard clip examples:
    >>> mismatch_pos([(5, 65), (0, 33)], '31A1')
    [31]
    >>> mismatch_pos([(0, 37), (5, 61)], '37')
    []
    '''
    mismatches = re.split('[^0-9, \^]+', mismatch_tag) #split by letters
    mis_pos = []
    index = 0    
    
    del_pos = 0 # keep track of pos # before deletion
    prev_del = False    
    
    for i in range(len(mismatches)-1):  
        # SNP in the first position         
        if mismatches[i] == '':
            mis_pos.append(0)
            continue
        # Ignore deletions, add pos num to subsequent mismatches
        if '^' in mismatches[i]:
            del_pos += int(mismatches[i][:-1])
            prev_del = True
            continue 
        
        mis_pos.append(int(mismatches[i]))
        
        # If prev pos contains deletion, add prev pos to current as we're ignoring deletions
        if prev_del:
            mis_pos[-1] += del_pos
            prev_del = False
            del_pos = 0

        if len(mismatches) >= 2:
            if i != 0 and len(mis_pos) > 1:
                # Need to add 1 to positions for correct indexing (except for first position)
                mis_pos[-1] += mis_pos[-2] + 1
    
    # Incorporate insertions 
    insert = 0
    for i in cigar:
        if i[0] == 1:
            # If there's multiple insertions
            for j in range(i[1]):
                mis_pos.append(insert)
                insert += 1
        elif i[0] == 0:
            # Keep track of positional num
            insert += i[1]
        else:
            pass
        
    mis_pos = list(set(mis_pos)) # Incase insertion and SNP share same position and there's repeat pos
    mis_pos.sort()
    
    return mis_pos


def query_seq_pos(cigar, readLength):
    '''(list of tuples, int) -> tuples
    Return tuple of seq position excluding soft clips and hard clips (0-based).
    
    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip
    
    >>> query_seq_pos([(4, 73), (0, 20), (4, 5)], 98)
    (73, 93)
    >>> query_seq_pos([(4, 6), (0, 92)], 98)
    (6, 98)
    >>> query_seq_pos([(0, 23), (4, 75)], 98)
    (0, 23)
    >>> query_seq_pos([(0, 37), (5, 61)], 98)
    (0, 37)
    '''
    start = 0 
    end = readLength
    
    if cigar[0][0] == 4 or cigar[0][0] == 5:
        start += cigar[0][1]
        
    if cigar[-1][0] == 4 or cigar[-1][0] == 5:
        end -= cigar[-1][1]

    return start, end


def consensus_maker(readList, readLength, cutoff):
    '''(list, int, int) -> str
    Return consensus sequence (without soft-clips) and quality score consensus.
    
    Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position. 
    
    - Add N's for soft clipped regions so it aligns with full length sequences
    
    - At each position, add quality score to list corresponding to nucleotide. 
      Take max quality score of nucleotide with highest frequency
    '''
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = array.array('B')
    
    mismatch_pos_lst = []
    
    for read in readList:
        mismatch_pos_lst.append(mismatch_pos(read.cigar, read.get_tag('MD')))
            
    for i in range(readLength):
        position_score = [0, 0 ,0, 0, 0] # A, C, G, T, N 
        quality_score = [[], [], [], [], []] 
        phred_fail = 0
        
        ### HOW MANY Ns ARE IN THE FINAL CONSENSUS AND HOW MANY TIE BREAKING EVENTS? ###
        
        for j in range(len(readList)):
            # Find position of sequence without soft clips
            query_pos = query_seq_pos(readList[j].cigar, readLength)
            # if seq length < or > region of query seq, add 1 to N and set qual score 0
            if i < query_pos[0] or i >= query_pos[1]:
                position_score[4] += 1
                quality_score[4].append(0)
                continue

            # Phred filter mismatch positions            
            if i in mismatch_pos_lst[j]:
                if readList[j].query_alignment_qualities[i] < 30: # Phred cutoff of 30
                    phred_fail += 1
                    continue
            
            i = i - query_pos[0] # account for seq with soft/hard clips
            # index subtract clipped bps to iterate through sequence
            # (e.g. 2S96M -> indexes 0 and 1 are N,
            # but at index 2 actual position in query seq is 0)

            # If pass filter, add 1 to nuc
            nuc = readList[j].query_alignment_sequence[i]
            nuc_index = nuc_lst.index(nuc)
        
            position_score[nuc_index] += 1  
            quality_score[nuc_index].append(readList[j].query_alignment_qualities[i])
            
            i = i + query_pos[0]                        


        try:
            # Find most common nuc #
            max_nuc_pos = [f for f, k in enumerate(position_score) if k == max(position_score)]
            # If there's more than one max, randomly select nuc
            max_nuc = max_nuc_pos[randint(0, len(max_nuc_pos)-1)]

	    # Median quality score (round to larger value if qual score is a decimal)
            max_qual = math.ceil(statistics.median(quality_score[max_nuc]))
        
            # frequency of nuc at position > cutoff 
            if max(position_score)/(len(readList) - phred_fail) > cutoff:
                consensus_read += nuc_lst[max_nuc]
                quality_consensus.append(max_qual)
            else:
                raise ValueError
                            
        except:
            # For cases when # matches fail phred > # reads
            consensus_read += 'N'
            quality_consensus.append(0)            

    return consensus_read, quality_consensus


def chr_arm_pos(chr_lst, chr_len):
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
    chr_arm_coor = collections.OrderedDict()
    
    if 'chrM' in chr_lst:
        chr_arm_coor['chrM'] = (0, chr_len[chr_lst.index('chrM')])
    
    filepath = os.path.abspath(inspect.getfile(inspect.currentframe())).rsplit('/', 1)[0]
    with open(filepath + '/cytoBand.txt') as f:
        next(f) # Skip header
        for line in f:
            chr_arm = line.split('\t')
            chr_key = '{}_{}'.format(chr_arm[0], chr_arm[3])
            start = int(chr_arm[1]) # python is 0-based (start is usually 1)
            end = int(chr_arm[2])
            chr_val = (start, end)

            chr_arm_coor[chr_key] = chr_val
    
    return chr_arm_coor


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--cutoff", action = "store", dest="cutoff", help="nucleotide base % cutoff", required = True)
    parser.add_argument("--Ncutoff", action = "store", dest="Ncutoff", help="N % cutoff", required = True)
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output SSCS BAM file", required = True)
    args = parser.parse_args()
    
    # need separate outfile argument as path is diff from input
    
    start_time = time.time()
    
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")    
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}_stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    
    # setup fastq files
    fastqFile1 = open('{}.r1.fastq'.format(args.outfile.split('.sscs')[0]), 'w')
    fastqFile2 = open('{}.r2.fastq'.format(args.outfile.split('.sscs')[0]), 'w')
    
    
    time_tracker = open('{}_time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    
    bam_dict = collections.OrderedDict() # dict subclass that remembers order entries were added
    tag_dict = collections.defaultdict(int)
       
    unmapped = 0
    unpaired = 0
    counter = 0  
    singletons = 0
    SSCS_reads = 0      
        
    chrm = [x['SN'] for x in bamfile.header['SQ']]
    chr_len = [x['LN'] for x in bamfile.header['SQ']]
    
    chr_arm_coor = chr_arm_pos(chrm, chr_len)
    #print(chr_arm_coor)
    
    for x in chr_arm_coor.keys():
        #print(x)
        bamLines = bamfile.fetch(reference = x.split('_')[0], start = chr_arm_coor[x][0], end =  chr_arm_coor[x][1]) # genomic start and end 0-based       
        # Create dictionary for each chrm
        for line in bamLines:
            counter += 1
            strand = 'fwd'
            if line.is_reverse:
                strand = 'rev'
                
            read = 'R1'
            if line.is_read2:
                read = 'R2'            
            
            tag = '{}_{}_{}_{}_{}_{}'.format(line.qname.split("|")[1], # mol barcode
                                          line.reference_name, # chr num
                                          line.reference_start, # start R1 (0-based)
                                          line.next_reference_start, # start R2
                                          strand, # strand direction
                                          read # read num
                                          )      
            
            if line.is_unmapped:
                unmapped += 1
                continue             
            elif not line.is_paired:
                unpaired += 1
                continue
            
            try:
                line.get_tag('MD')
            except:
                print(line)
                continue
                
            tag_dict[tag] += 1     
            
            if tag not in bam_dict:
                bam_dict[tag] =[line]
        
            else:
                bam_dict[tag].append(line)             
                
                
        # ===== Create consenus seq for reads in each chrm arm and reset =====
        
        tag_keys = bam_dict.keys()            
        for i in tag_keys:
            if tag_dict[i] < 2:
                singletons += 1
                singleton_bam.write(bam_dict[i][0])
            else:
                readLength = max(collections.Counter(i.query_alignment_length for i in bam_dict[i]))

                SSCS = consensus_maker(bam_dict[i], readLength, float(args.cutoff))
                
                ## NEED TO FIX FOR ALL CONSENSUS MAKING!!!!! Need to divide by total number of bases                 
                if SSCS[0].count('N')/len(SSCS[0]) > float(args.Ncutoff):
                    continue
		
                # write as bam
                SSCS_read = create_aligned_segment(bam_dict[i], SSCS[0], SSCS[1])
                SSCS_bam.write(SSCS_read)

		# write as fastq file
                if "R1" in i:
		    # fastq format: query_name, consensus_seq, qual score
                    fastqFile1.write('@:{}\n{}\n+\n{}\n'.format(SSCS_read.qname, SSCS[0], pysam.qualities_to_qualitystring(SSCS[1])))
                else:
                    fastqFile2.write('@:{}\n{}\n+\n{}\n'.format(SSCS_read.qname, SSCS[0], pysam.qualities_to_qualitystring(SSCS[1])))
                
                
                SSCS_reads += 1
	   
        # reset dictionary            
        bam_dict = collections.OrderedDict() # dict subclass that remembers order entries were added        
        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')
            
        except:
            continue
    
    # ===== write tag family size dictionary to file ===== 
    # (key = tags, value = int [number of reads in that family]) 
    import pickle
    tag_file = open(args.outfile.split('.sscs')[0] + '.read_families.txt', 'ab+')
    pickle.dump(tag_dict, tag_file)
    tag_file.close()
    
    summary_stats = '''Total reads: {} \n
Unmapped reads: {} \n
Unpaired reads: {} \n
SSCS reads: {} \n
Singletons: {} \n
'''.format(counter, unmapped, unpaired, SSCS_reads, singletons)

    stats.write(summary_stats)
    
    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    singleton_bam.close()
    
    
    # ===== Create tag family size plot =====
    # Count number of families containing each number of read (e.g. Counter({1: 3737, 32: 660... -> 3737 families are singletons)
    fam_per_read_group = collections.Counter([i for i in tag_dict.values()])
    lst_fam_per_read = list(fam_per_read_group.items()) # conver to list
    
    total_reads = sum(tag_dict.values())
    # Multiply number of families by read num to get total number of reads in that read group, divide it by total reads to obtain read fraction
    read_fraction = [(i*j)/total_reads for i,j in lst_fam_per_read] 
    
    plt.bar(list(fam_per_read_group), read_fraction)
    #plt.locator_params(axis = 'x', nbins=lst_fam_per_read[-1][0]//5)
    plt.xlabel('Tag family size (# of reads per family)')
    plt.ylabel('Fraction of total reads')
    
    plt.savefig(args.outfile.split('.sscs')[0]+'_tag_fam_size.png')      
    

###############################
##           Main            ##
###############################
if __name__ == "__main__": 
    import time
    start_time = time.time()
    main()  
    print((time.time() - start_time)/60)  


# use python, write to temp files, then merge together using samtools or bamtools and then delete temp files