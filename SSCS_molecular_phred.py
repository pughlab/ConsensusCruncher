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
def mismatch_pos(cigar, mismatch_tag):
    '''(list, str) -> lst
    Return 0-based index of mismatch positions in sequence (including insertions,
    ignoring deletions). List of positions to check phred quality score.
    
    cigar: soft clips, hard clips, insertions, deletions
    mismatch_tag: position of mismatch
    
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
    #mismatch_tag = mismatch_tag.split('N0')[0]
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
    
    ## NEED TO RE-WRITE - have to offset mismatch pos if there's insertion!!!!    
    
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
    
    Soft clips are shown in seq, hard clips are not.
    
    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip
    
    >>> query_seq_pos([(4, 73), (0, 20), (4, 5)], 98)
    (73, 93)
    >>> query_seq_pos([(5, 78), (0, 15), (4, 5)], 98)
    (0, 15)
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
    '''(list, int, int) -> strq
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
        #quality_score = [[], [], [], [], []] 
        #quality_score = [0, 0, 0, 0, 0]
        phred_fail = 0
        
        ### HOW MANY Ns ARE IN THE FINAL CONSENSUS AND HOW MANY TIE BREAKING EVENTS? ###
        
        for j in range(len(readList)):
            # === Find position of sequence without soft clips ===
            query_pos = query_seq_pos(readList[j].cigartuples, readLength)
            # if seq length < or > region of query seq, add 1 to N and set qual score 0
            if i < query_pos[0] or i >= query_pos[1]:
                position_score[4] += 1
                #quality_score[4].append(0)
                continue

            # === Phred filter mismatch positions ===     
            if i in mismatch_pos_lst[j]:
                if readList[j].query_alignment_qualities[i] < 30: # Phred cutoff of 30
                    phred_fail += 1
                    continue
            
            i = i - query_pos[0] # account for seq with soft/hard clips
            # index subtract clipped bps to iterate through sequence
            # (e.g. 2S96M -> indexes 0 and 1 are N,
            # but at index 2 actual position in query seq is 0)
            # IF you use query_sequence (which includes soft clips), then you don't need to offset the query pos => however, it will throw off consensus making between those with and without soft clips (e.g. readA: ACGTT, readB: TGACGTT (TG soft clips) but pos based comparison will show they're off) 

            # If pass filter, add 1 to nuc
            nuc = readList[j].query_alignment_sequence[i]
            nuc_index = nuc_lst.index(nuc)
        
            position_score[nuc_index] += 1  
            #quality_score[nuc_index] += 1
            #quality_score[nuc_index].append(readList[j].query_alignment_qualities[i])
            
            i = i + query_pos[0]                        


        try:
            # Find most common nuc #
            max_nuc_pos = [f for f, k in enumerate(position_score) if k == max(position_score)]
            # If there's more than one max, randomly select nuc
            max_nuc_pos = max_nuc_pos[randint(0, len(max_nuc_pos)-1)]
            
            # === Molecular Phred Quality Score ===
            # error = num variant bases (not most freq base)
            error_bases = position_score[:max_nuc_pos] + position_score[(max_nuc_pos +1):]
            # probability of observed bases
            P = sum(error_bases)/sum(position_score)
            if P == 0:
                Q = 62
            else:
                Q = round(-10 * math.log10(P))
                
            consensus_read += nuc_lst[max_nuc_pos]
            quality_consensus.append(Q)
        
            # frequency of nuc at position > cutoff 
            #if max(position_score)/(len(readList) - phred_fail) > cutoff:
                #consensus_read += nuc_lst[max_nuc_pos]
                #quality_consensus.append(max_qual)
            #else:
                #raise ValueError
                            
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


def reverse_seq(seq):
    '''(str) -> str
    Return reverse compliment of sequence (used for writing rev comp sequences to fastq files).
    
    >>> reverse_seq('TCAGCATAATT')
    'AATTATGCTGA'
    >>> reverse_seq('ACTGNN')
    'NNCAGT'
    '''
    rev_comp = ''
    nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in seq:
        rev_comp = nuc[base] + rev_comp
        
    return rev_comp


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--cutoff", action = "store", dest="cutoff", help="nucleotide base % cutoff", required = True)
    parser.add_argument("--Ncutoff", action = "store", dest="Ncutoff", help="N % cutoff", required = True)
    parser.add_argument("--infile", action = "store", dest="infile", help="input BAM file", required = True)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="output SSCS BAM file", required = True)
    args = parser.parse_args()
    
    
    start_time = time.time()
    
    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)    
    doubleton_bam = pysam.AlignmentFile('{}.doubleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    
    # setup fastq files
    fastqFile1 = open('{}.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    fastqFile2 = open('{}.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    
    doubleton_fastqFile1 = open('{}.doubleton.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    doubleton_fastqFile2 = open('{}.doubleton.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')    
    
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    
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
    
    chr_arm_coor = chr_arm_pos(chrm, chr_len)
    
    for x in chr_arm_coor.keys():
        # Create dictionary for each chrm
        chr_data = read_bam(bamfile, read_chr = x.split('_')[0], read_start = chr_arm_coor[x][0], read_end = chr_arm_coor[x][1]) # genomic start and end 0-based
        
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
        
        
        # ===== Create consenus seq for reads in each chrm arm and reset =====
        if bool(bam_dict):
            rand_key = choice(list(bam_dict.keys()))
            readLength = bam_dict[rand_key][0].infer_query_length()        
        
        written_pairs = []
        
        for readPair in trans_pair.keys():
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
                        aligned_seq = SSCS_read.query_alignment_sequence
                        if aligned_seq.count('N')/len(aligned_seq) > float(args.Ncutoff):
                            print('uh oh too many Ns')
                            print(aligned_seq)
                            print(SSCS_read)
                            continue                        

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
                        
                        # Use aligned sequence in case of soft clips
                        aligned_seq = SSCS_read.query_alignment_sequence
                        if aligned_seq.count('N')/len(aligned_seq) > float(args.Ncutoff):
                            print('uh oh too many Ns')
                            print(aligned_seq)
                            print(SSCS_read)
                            continue                        

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
    
    
    ## READ IN BAM FILE LATER AND CREATE THESE PLOTS
    # ===== Create quality score plot =====
    
    # Should we take average quality score of every read??????
    import itertools
    
    qual_dist_lst = [list(i[0]) for i in quality_dict.values()]
    qual_dist = list(itertools.chain.from_iterable(qual_dist_lst))
    count_qual = collections.Counter(qual_dist).most_common()
    count_qual_sorted = sorted(count_qual, key=lambda tup: tup[0])
    x = [x for x,y in count_qual_sorted]
    y = [y for x,y in count_qual_sorted]    
    
    fig = plt.figure(figsize=(7.195, 3.841), dpi=100)  
    ax = fig.add_subplot(111)
    
    plt.plot(x, y, "-o", markersize=np.sqrt(5), linewidth = 0.5)
    plt.yscale('log')
    plt.xticks(np.arange(0, max(x)+5, 5.0), fontsize = 6)
    plt.yticks(fontsize = 6)    
    plt.xlabel('Phred Quality Score', fontsize = 8)
    plt.ylabel('Number of bases', fontsize = 8)
    
    for xy in zip(x, y):                                     
        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data', fontsize = 4)     
    plt.grid()     
    
    plt.savefig(args.outfile.split('.sscs')[0]+'.phred_density.png', dpi = 1000) 
    
    
    # ===== Create quality score plot for each tag family =====
    
    # Figure out how to plot tag_quality_dict....
    
    # PLOT BASES -> density plot phred vs bases
    
    x = [[i]*len(tag_quality_dict[i]) for i in tag_quality_dict]
    y = list(tag_quality_dict.values())
    
    x = list(itertools.chain.from_iterable(x))
    y = list(itertools.chain.from_iterable(y))
    
    plt.figure(figsize=(7.195, 3.841), dpi=100)
    
    plt.plot(x, y, "o", markersize=np.sqrt(5))
    
    plt.xticks(np.arange(0, max(x)+1, 5.0), fontsize = 6)
    plt.yticks(np.arange(0, max(y)+5, 5.0), fontsize = 6)
    
    plt.ylabel('Average Molecular Q / Read', fontsize = 8)
    plt.xlabel('Tag family size', fontsize = 8)    
    
    #tick.label.set_fontsize(8)
    
    plt.grid() 
    
    plt.savefig(args.outfile.split('.sscs')[0]+'.tag_fam_quality.png', dpi = 1000)
    

###############################
##           Main            ##
###############################
if __name__ == "__main__": 
    import time
    start_time = time.time()
    main()  
    print((time.time() - start_time)/60)  


# use python, write to temp files, then merge together using samtools or bamtools and then delete temp files




#if line.is_unmapped:
    #unmapped += 1
    #mate_read = 'R1'
    #if read == 'R1':
        #mate_read = 'R2'
    #unmapped_key += [tag[:-2] + mate_read]
    
    ##try:
        ##unmapped_reads += [bamfile.mate(line)]
    ##except:
        ##print('uh oh, this read has no mate! {}'.format(line))
    ##print(line)
    ##print(tag)

    ##print(bamfile.mate(line))
    ##break 
    #continue
#elif tag in unmapped_key:
    #pair_unmapped += 1
    #continue
#elif line.reference_start == line.next_reference_start:
    #try:
        #if bamfile.mate(line).is_unmapped:
            #pair_unmapped += 1
            #continue
    #except:
        #print(line)
        #pair_unmapped += 1
        #continue
        
    ##print(line)
    ###print(bamfile.mate(line))
    ##if bamfile.mate(line).is_unmapped:
        ##continue
    ##else:
        ##print('uh oh')



        
        ##bamfile = pysam.AlignmentFile('/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/subsampled_bam/epic_md_tag/SWID_4627086_EPICT_057_nn_C_PE_202_TS_NoGroup_160621_D00331_0196_AC900FANXX_CTAAGTGG_L007.processed.bam', 'rb')
        #bamfile = pysam.AlignmentFile('/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test5/MEM-001-KRAS.bam', 'rb')
        
        #a = read_bam(bamfile)
        
        #readList = a[0]['TTCC_12_25398453_12_25398323_rev_R2']
        
        #for i in readList:
            #print(i)
            
        
        #if bool(a[0]):
            #rand_key = choice(list(a[0].keys()))
            #readLength = a[0][rand_key][0].infer_query_length()     
        
        #readList = readList[14:16]
        #print(consensus_maker(readList, readLength, 0.7))
        
        
        
        
        # =============
        
        # consensus_maker fx
        
        
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