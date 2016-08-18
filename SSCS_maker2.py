#!/usr/bin/env python3

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
from random import *
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import statistics
import math
import time
import inspect
import os

from consensus_helper import *

###############################
##         Functions         ##
###############################

def mismatch_pos(cigar, mismatch_tag, readLength):
    '''(list, str) -> lst
    Return 0-based index of mismatch positions in sequence (including insertions,
    ignoring deletions) - list of variant positions.
    
    === Background Info ===
    == CIGAR ==  
    - soft clips, hard clips, insertions, deletions
    - Soft clips are shown in seq, hard clips are not.
    
    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip    
    
    == Mismatch_tag/MD ==
    - position of mismatch 
    ** NOTE: MD only represent information about read aligned to reference, don't take into consideration insertions
    
    e.g. mismatch_pos([(0, 28), (1, 1), (0, 69)], '19T74G2')
    total bases: cigar = 98, MD = 97
    
    ========================
    ** Mismatch positions need to account for both cigar and MD (offset MD positions when insertions are present)
    
    'N0N0N0' seq due to non-overlap region of read with ref (end of ref seq)
    
    === Test cases ===
    # 14S29M2D30M50S, cigar = 59, MD = 59 
    >>> mismatch_pos([(4, 14), (0, 29), (2, 2), (0, 30), (4, 50)], '13c0g14^ga3g23N0N0N0', 123)
    [14, 27, 28, 46, 70, 71, 72, 73]
    
    ((14, 73), [13, 14, 32, 56, 57, 58])
    
    # 47M3D27M6I13M30S
    >>>mismatch_pos([(0, 47), (2, 3), (0, 27), (1, 6), (0, 13), (4, 30)], '19a27^ggc12g25N0N0', 123)
    [0, 19, 59, 74, 75, 76, 77, 78, 79, 91, 92, 93]
    
    ((0, 93), [19, 59, 74, 75, 76, 77, 78, 79, 91, 92])
    
    # 35M4I16M5I32M31H <- realistically, hard clips are filtered out
    >>>mismatch_pos([(0, 35), (1, 4), (0, 16), (1, 5), (0, 32), (5, 31)], '3g10g10g1g15c38N0', 92)
    [0, 3, 14, 25, 27, 35, 36, 37, 38, 47, 55, 56, 57, 58, 59, 91, 92]
    
    ((0, 92), [3, 14, 25, 27, 35, 36, 37, 38, 47, 55, 56, 57, 58, 59, 60, 91])
    
    8S115M
    >>>mismatch_pos([(4, 8), (0, 115)], '0N0N113', 123)
    [8, 8, 9, 123]
    
    E.g. mismatch_pos([(0, 19), (2, 1), (0, 79)], '19^A8G70', 98)
    [27]
    ^ = deletion 
    19 matches, followed by 1bp deletion, 8 matches, with a G on the reference which is different form aligned read base, and 70 matches
    
    First mismatch "G' at position 19 + 8
    '''
    # Adjust for softclips - CHECK TO SEE IF INDEXED CORRECTLY
    start = 0
    end = readLength
    
    if cigar[0][0] == 4:
        start += cigar[0][1]
    if cigar[-1][0] == 4:
        end -= cigar[-1][1]    
    
    # Identify insertions in cigar tuples [(operation, value)]  
    # Insertion index corresponds to ALIGNED seq (soft clips excluded)
    inserts = []
    index = start
    for operation in cigar:
        if operation[0] == 1:    
            inserts.append((index, operation[1]))
            index += operation[1]
        elif operation[0] == 0:
            index += operation[1]
        else:
            pass    
        
    # Mismatch position for aligned sequence (excluding soft clips)
    mismatches = re.split('[^0-9, \^]+', mismatch_tag) #split by letters
    mis_pos = [start, end]
    index = start  
    
    inserts = sorted(inserts)
    old_inserts = []
    
    for base in range(len(mismatches)-1):
        # SNP in first pos
        if mismatches[base] == '':
            mis_pos.append(0)
        elif '^' in mismatches[base]:
            index += int(mismatches[base][:-1])
        else:
            index += int(mismatches[base])
            mis_pos.append(index)
            index += 1
            for i in inserts:
                if index > i[0] and i not in old_inserts:
                    mis_pos[-1] = mis_pos[-1] + i[1]
                    for val in range(i[1]):
                        mis_pos.append(i[0] + val)
                        index += 1
                    old_inserts.append(i)
    
    return sorted(mis_pos)


def consensus_maker(readList, readLength, cutoff):
    '''(list, int, int) -> strq
    Return consensus sequence (without soft-clips) and quality score consensus.
    
    Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position. 
    
    - Add N's for soft clipped regions so it aligns with full length sequences
    
    - At each position, add quality score to list corresponding to nucleotide. 
      Take max quality score of nucleotide with highest frequency
      
    Bases below the phred quality cutoff (0.7) are excluded from consenus making 
    '''
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = array.array('B')
    
    mismatch_pos_lst = []
    
    for read in readList:
        mismatch_pos_lst.append(mismatch_pos(read.cigartuples, read.get_tag('MD'), read.infer_query_length()))

    for i in range(readLength):
        position_score = [0, 0 ,0, 0, 0] # A, C, G, T, N
        quality_score = [[], [], [], [], []] 
        phred_fail = 0
        
        for j in range(len(readList)):
            # === Add N and set qual score to 0 for soft clipped regions ===
            if i < mismatch_pos_lst[j][0] or i >= mismatch_pos_lst[j][-1]:
                position_score[4] += 1
                quality_score[4].append(0)
            else:
                # === Phred filter mismatch positions - Phred cutoff of 30 ===
                # need to either offset i or add start to the mismatch positions
                if i in mismatch_pos_lst[j][1:-1]:
                    if readList[j].query_qualities[i] < 30: 
                        phred_fail += 1
                        continue
                
                nuc = readList[j].query_sequence[i]
                nuc_index = nuc_lst.index(nuc)
                position_score[nuc_index] += 1
                quality_score[nuc_index].append(readList[j].query_qualities[i])
                
        # Find most common nuc 
        # there's scenarios were phred fail == len(readList) resulting in empty quality score list, this also causes zero division error 
        try:
            max_nuc_index = [f for f, k in enumerate(position_score) if k == max(position_score)]
            # If there's more than one max, randomly select nuc
            max_nuc = max_nuc_index[randint(0, len(max_nuc_index)-1)]
            # Median quality score (round to larger value if qual score is a decimal)
            import statistics
            import math
            med_qual = math.ceil(statistics.median(quality_score[max_nuc]))

            # frequency of nuc at position > cutoff
            if position_score[max_nuc]/(len(readList) - phred_fail) > cutoff:
                consensus_read += nuc_lst[max_nuc]
                quality_consensus.append(med_qual)
            else:
                raise ValueError                
    
        except ValueError:
            consensus_read += 'N'
            quality_consensus.append(0)        
        
        #max_nuc_index = [f for f, k in enumerate(position_score) if k == max(position_score)]
        ## If there's more than one max, randomly select nuc
        #max_nuc = max_nuc_index[randint(0, len(max_nuc_index)-1)]
        ## Median quality score (round to larger value if qual score is a decimal)
        #import statistics
        #import math
        #try:
            #med_qual = math.ceil(statistics.median(quality_score[max_nuc]))
        #except:
            #print(position_score)
            #print(quality_score)
            #print(quality_score[max_nuc])
            #print(phred_fail)
            #print(i)
            #print(readList[0])
            #print(readList[1])
            #return 'hi'
        
        ## frequency of nuc at position > cutoff
        ##if len(readList) != phred_fail:
        #try:
            ##print(position_score[max_nuc])
            ##print((len(readList) - phred_fail))
            #if position_score[max_nuc]/(len(readList) - phred_fail) > cutoff:
                #consensus_read += nuc_lst[max_nuc]
                #quality_consensus.append(med_qual)
            #else:
                #raise ValueError
        #except ValueError:
            #consensus_read += 'N'
            #quality_consensus.append(0)
        #print(i)
        #print(consensus_read)
        
    return consensus_read, quality_consensus 


def chr_arm_pos(chr_lst, chr_len, bedfile = None):
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
    
    if bedfile == None:
        filepath = os.path.abspath(inspect.getfile(inspect.currentframe())).rsplit('/', 1)[0]
        bedfile = filepath + '/cytoBand.txt'
    
    with open(bedfile) as f:
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
    #parser.add_argument("--Ncutoff", action = "store", dest="Ncutoff", help="N % cutoff", required = True)
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output SSCS BAM file", required = True)
    parser.add_argument("--bedfile", action = "store", dest="bedfile", help="input bedfile", required = False)
    args = parser.parse_args()
    
    start_time = time.time()
    
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)    
    doubleton_bam = pysam.AlignmentFile('{}.doubleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    
    # set up fastq files
    fastqFile1 = open('{}.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    fastqFile2 = open('{}.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    
    doubleton_fastqFile1 = open('{}.doubleton.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    doubleton_fastqFile2 = open('{}.doubleton.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')    
    
    # set up time tracker
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    
    # ===== Initialize dictionaries =====
    trans_pair = collections.OrderedDict()
    trans_dict = collections.OrderedDict()
    
    tag_quality_dict = collections.defaultdict(list)
    quality_dict = collections.defaultdict(list)
    
    unmapped = 0
    unmapped_flag = 0
    bad_reads = 0 # secondary/supplementary reads
    counter = 0  
    doubletons = 0
    singletons = 0
    SSCS_reads = 0    
    
    chrm = [x['SN'] for x in bamfile.header['SQ']]
    chr_len = [x['LN'] for x in bamfile.header['SQ']]
    
    if 'args.bedfile' in locals():
        chr_arm_coor = chr_arm_pos(chrm, chr_len, args.bedfile)
    else:   
        chr_arm_coor = chr_arm_pos(chrm, chr_len)
    
    for x in chr_arm_coor.keys():
        # Create dictionary for each chrm
        #print(trans_pair)        
        ## carry over information from trans dict from one chr to another
        chr_data = read_bam(bamfile, read_chr = x.split('_')[0], read_start = chr_arm_coor[x][0], read_end = chr_arm_coor[x][1], trans_pair = trans_pair, trans_dict = trans_dict) # genomic start and end 0-based
        
        # Dictionary is reset each loop to contain only data from given chr coor
        bam_dict = chr_data[0]
        tag_dict = chr_data[1]
        paired_dict = chr_data[2]
        
        trans_dict.update(chr_data[3])
        trans_pair.update(chr_data[4])
        
        counter += chr_data[5]
        unmapped += chr_data[6]
        unmapped_flag += chr_data[7]
        bad_reads += chr_data[8]    
        
        #print(paired_dict['CTTC_1_116099948_1_116100023_neg_83_163'])
        #print(tag_dict['CTTC_1_116100023_1_116099948_rev_R1'])
        
    
        # ===== Create consenus seq for reads in each chrm arm and reset =====
        if bool(bam_dict):
            rand_key = choice(list(bam_dict.keys()))
            readLength = bam_dict[rand_key][0].infer_query_length()        
        
        written_pairs = []
        
        for readPair in trans_pair.keys():
            #if readPair == 'ACGC_12_120699897_12_120699962_neg_83_163':
                #print(trans_pair[readPair])
            if len(trans_pair[readPair]) == 2:
                for tag in trans_pair[readPair]:
                    # Check for singletons
                    if tag_dict[tag] < 2:
                        singletons += 1
                        singleton_bam.write(trans_dict[tag][0])
                    else:
                        SSCS = consensus_maker(trans_dict[tag], readLength, float(args.cutoff))
                        
                        quality_dict[tag] += [SSCS[1]]
                        tag_quality_dict[tag_dict[tag]] += [round(np.mean(SSCS[1]))]
                        
                        SSCS_read = create_aligned_segment(trans_dict[tag], SSCS[0], SSCS[1])
                        
                        #new_tag = sscs_qname(tag, SSCS_read.flag)
                        query_name = readPair + ':' + str(tag_dict[tag])   
                        SSCS_read.query_name = query_name
                        
                        # Use aligned sequence in case of soft clips
                        #aligned_seq = SSCS_read.query_alignment_sequence
                        #if aligned_seq.count('N')/len(aligned_seq) > float(args.Ncutoff):
                            #print('uh oh too many Ns')
                            #print(aligned_seq)
                            #print(SSCS_read)
                            #continue                        

                        if tag_dict[tag] == 2:
                            doubletons += 1
                            doubleton_bam.write(SSCS_read)
                        else:
                            SSCS_bam.write(SSCS_read)    
                            SSCS_reads += 1
                
                        # ===== write as fastq file =====                        
                        #if SSCS_read.is_reverse and SSCS_read.is_read1:
                        #fastq_seq = SSCS_read.query_sequence.decode("utf-8") 
                        fastq_seq = SSCS_read.query_sequence 
                        
                        
                        if 'rev' in tag:
                            fastq_seq = reverse_seq(fastq_seq)
                            fastq_qual = pysam.qualities_to_qualitystring(reversed(SSCS_read.query_qualities))
                        else:
                            fastq_qual = pysam.qualities_to_qualitystring(SSCS_read.query_qualities)
                        
                        if tag_dict[tag] == 2:
                            if 'R1' in tag:
                                doubleton_fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                            else:
                                doubleton_fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))                            
                        else:
                            if 'R1' in tag:
                                fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                            else:
                                fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                written_pairs.append(readPair)
                
        for read in written_pairs:
            trans_pair.pop(read)        
        
        for readPair in paired_dict.keys():          
            # Check pairing
            if len(paired_dict[readPair]) != 2:
                print('uh oh pairing problem!!!')
                print(chr_arm_coor[x][0])
                print(chr_arm_coor[x][1])
                print(readPair)
                print(paired_dict[readPair])
                print(bam_dict[paired_dict[readPair][0]][0])
                m = bamfile.mate(bam_dict[paired_dict[readPair][0]][0])
                print(m)
                
                strand = 'fwd'
                if m.is_reverse:
                    strand = 'rev'
                    
                read = 'R1'
                if m.is_read2:
                    read = 'R2'            
                mtag = '{}_{}_{}_{}_{}_{}_{}'.format(m.qname.split("|")[1], # mol barcode
                                                      m.reference_id, # chr num
                                                      m.reference_start, # start R1 (0-based)
                                                      m.next_reference_id,
                                                      m.next_reference_start, # start R2
                                                      strand, # strand direction
                                                      read # read num
                                                      )             
                print(mtag)
                print(sscs_qname(mtag, m.flag))                
                return readPair                
            else:
                for tag in paired_dict[readPair]:
                    # Check for singletons
                    if tag_dict[tag] < 2:
                        singletons += 1
                        singleton_bam.write(bam_dict[tag][0])
                    else:
                        SSCS = consensus_maker(bam_dict[tag], readLength, float(args.cutoff))
                        
                        quality_dict[tag] += [SSCS[1]]
                        tag_quality_dict[tag_dict[tag]] += [round(np.mean(SSCS[1]))]
                        
                        SSCS_read = create_aligned_segment(bam_dict[tag], SSCS[0], SSCS[1])
                        
                        #new_tag = sscs_qname(tag, SSCS_read.flag)
                        query_name = readPair + ':' + str(tag_dict[tag])   
                        SSCS_read.query_name = query_name
                        
                        ## Use aligned sequence in case of soft clips
                        #aligned_seq = SSCS_read.query_alignment_sequence
                        #if aligned_seq.count('N')/len(aligned_seq) > float(args.Ncutoff):
                            #print('uh oh too many Ns')
                            #print(aligned_seq)
                            #print(SSCS_read)
                            #continue                        

                        if tag_dict[tag] == 2:
                            doubletons += 1
                            doubleton_bam.write(SSCS_read)
                        else:
                            SSCS_bam.write(SSCS_read)    
                            SSCS_reads += 1
                
                        # ===== write as fastq file =====                        
                        #if SSCS_read.is_reverse and SSCS_read.is_read1:
                        #fastq_seq = SSCS_read.query_sequence.decode("utf-8") 
                        fastq_seq = SSCS_read.query_sequence 
                        
                        
                        if 'rev' in tag:
                            fastq_seq = reverse_seq(fastq_seq)
                            fastq_qual = pysam.qualities_to_qualitystring(reversed(SSCS_read.query_qualities))
                        else:
                            fastq_qual = pysam.qualities_to_qualitystring(SSCS_read.query_qualities)
                        
                        if tag_dict[tag] == 2:
                            if 'R1' in tag:
                                doubleton_fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                            else:
                                doubleton_fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))                            
                        else:
                            if 'R1' in tag:
                                fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                            else:
                                fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                        
        import pickle
        read_dict_file = open(args.outfile.split('.sscs')[0] + '.read_dict.txt', 'ab+')
        pickle.dump(bam_dict, read_dict_file)
        read_dict_file.close()
        
        # reset dictionary            
        bam_dict = collections.OrderedDict() # dict subclass that remembers order entries were added        
        paired_dict = collections.OrderedDict()
        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')
            
        except:
            continue
    
    # ===== write tag family size dictionary to file ===== 
    # (key = tags, value = int [number of reads in that family]) 
    tag_file = open(args.outfile.split('.sscs')[0] + '.read_families.txt', 'ab+')
    pickle.dump(tag_dict, tag_file)
    tag_file.close()
    
    summary_stats = '''Total reads: {} \n
Unmapped reads: {} \n
Unmapped flag reads: {} \n
Secondary/Supplementary reads: {} \n
SSCS reads: {} \n
singletons: {} \n
doubletons: {} \n
'''.format(counter, unmapped, unmapped_flag, bad_reads, SSCS_reads, singletons, doubletons)

    stats.write(summary_stats)
    
    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    doubleton_bam.close()
    fastqFile1.close()
    fastqFile2.close()
    doubleton_fastqFile1.close()
    doubleton_fastqFile2.close()      
    
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