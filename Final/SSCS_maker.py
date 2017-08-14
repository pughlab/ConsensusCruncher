#!/usr/bin/env python3

###############################################################
#
#                           SSCS Maker
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
#
# Concepts:
#    - Read family: reads that share the same molecular barcode, genome
#                   coordinates for Read1 and Read2, cigar string, strand, flag, and read number
#    - Singleton: a read family containing only one member (a single read)
#
# Usage: python3 SSCS_maker.py [--cutoff CUTOFF] [--infile INFILE] [--outfile OUTFILE]
#
# Arguments:
# --infile INFILE     input BAM file
# --outfile OUTFILE   output BAM file
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
###############################################################

import pysam # Need to install
import collections
import re
import array
from random import *
from itertools import chain
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import statistics
import math
import time

from consensus_helper import *

###############################
#          Functions          #
###############################


def consensus_maker(readList, cutoff):
    '''(list, int, int) -> str
    Return consensus sequence and quality score.

    Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position.
    - At each position, reads supporting each nucleotide is recorded along with the quality score corresponding to each
    nucleotide
    - Bases below the phred quality cutoff (Q30) are excluded from consensus making
    - The most frequent base is added to the consensus sequence, given that the proportion of reads supporting this base
    is greater than the cutoff
    - A molecular phred quality score (consensus quality score) is determined by taking the product of errors of the most
    frequent base
    - If a majority can't be determined (i.e. a tie with 2 maximums)
    '''
    # === Determine read length ===
    cigar_mode = statistics.mode([r.cigarstring for r in readList])
    cigar_mode_read = [r for r in readList if r.cigarstring == cigar_mode][0]
    readLength = cigar_mode_read.infer_query_length()

    # === Initialize counters ===
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = []
    proportion_scores = []

    for i in range(readLength):
        nuc_count = [0, 0, 0, 0, 0]  # A, C, G, T, N
        quality_score = [[], [], [], [], []]
        phred_fail = 0

        for j in range(len(readList)):
            # === Phred filter Q30 cut-off ===
            if readList[j].query_qualities[i] < 30:
                phred_fail += 1
            else:
                nuc = readList[j].query_sequence[i]
                nuc_index = nuc_lst.index(nuc)
                nuc_count[nuc_index] += 1
                quality_score[nuc_index].append(readList[j].query_qualities[i])

        # Find most frequent nuc (don't worry about ties (2 maxes) as it won't pass proportion cut-off and N will
        # be assigned)
        max_nuc = nuc_count.index(max(nuc_count))

        # === Molecular phred quality (consensus quality score) ===
        max_nuc_quality = quality_score[max_nuc]

        base_fail = False

        if max_nuc_quality != []:
            P = 1
            for Q in max_nuc_quality:
                P *= 10**(-(Q/10))

            mol_qual = round(-10 * math.log10(P))

            if mol_qual > 62:
                mol_qual = 62

        else:
            mol_qual = 0
            base_fail = True

        # === Check proportion cutoff ===
        # only make consensus if proportion of nuc is > cutoff (e.g. 70%) of reads
        phred_pass_reads = len(readList) - phred_fail
        if phred_pass_reads != 0:
            prop_score = nuc_count[max_nuc]/phred_pass_reads
            if prop_score >= cutoff:
                if base_fail == True:
                    # test to see if a position that has failed quality score making can still pass filters
                    print('base fail == True!!')
                consensus_read += nuc_lst[max_nuc]
                quality_consensus.append(mol_qual)
                proportion_scores.append(prop_score)
            else:
                base_fail = True
        else:
            base_fail = True

        # === Set base to N if no consensus could be made ===
        if base_fail:
            consensus_read += 'N'
            quality_consensus.append(mol_qual)
            proportion_scores.append(0)

    return consensus_read, quality_consensus, proportion_scores


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--cutoff", action="store", dest="cutoff", help="nucleotide base % cutoff", required=True)
    # parser.add_argument("--Ncutoff", action = "store", dest="Ncutoff", help="N % cutoff", required = True)
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output SSCS BAM file", required=True)
    parser.add_argument("--bedfile", action="store", dest="bedfile", help="input bedfile for data division", required=False)
    args = parser.parse_args()

    start_time = time.time()

    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    tag_bam = pysam.AlignmentFile('{}.tag_reads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    # set up time tracker
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')

    # ===== Initialize dictionaries =====
    read_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    pair_dict = collections.defaultdict(list)
    csn_pair_dict = collections.defaultdict(list)

    # ===== Initialize counters =====
    unmapped = 0
    unmapped_mate = 0
    multiple_mapping = 0  # secondary/supplementary reads
    counter = 0
    singletons = 0
    SSCS_reads = 0

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

        chr_data = read_bam(bamfile,
                            pair_dict=pair_dict,
                            read_dict=read_dict,
                            csn_pair_dict=csn_pair_dict,
                            tag_dict=tag_dict,
                            badRead_bam=badRead_bam,
                            duplex=None,
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

        # ===== Create consensus sequences as paired reads =====
        for readPair in list(csn_pair_dict.keys()):
            if len(csn_pair_dict[readPair]) == 2:
                for tag in csn_pair_dict[readPair]:
                    # === Check for singletons ===
                    if tag_dict[tag] == 1:
                        singletons += 1
                        # print(read_dict[tag])
                        singleton_bam.write(read_dict[tag][0])
                    else:
                        SSCS = consensus_maker(read_dict[tag], float(args.cutoff))

                        query_name = readPair + ':' + str(tag_dict[tag])
                        SSCS_read = create_aligned_segment(read_dict[tag], SSCS[0], SSCS[1], query_name)

                        # ===== Write consensus bam =====
                        SSCS_bam.write(SSCS_read)
                        SSCS_reads += 1

                    # Remove read from dictionary after writing
                    del read_dict[tag]

                # Remove key from dictionary after writing
                del csn_pair_dict[readPair]

        print(x)
        print(str((time.time() - start_time)/60))
        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')

        except:
            continue

    # === Check to see if there's remaining reads ===
    print('=== pair_dict remaining ===')
    if bool(pair_dict):
        for i in pair_dict:
            print(i)
            print('read remaining:')
            print(pair_dict[i][0])
            print('mate:')
            print(bamfile.mate(pair_dict[i][0]))
    print('=== read_dict remaining ===')
    if bool(read_dict):
        for i in read_dict:
            print(i)
            print('read remaining:')
            print(read_dict[i][0])
            print('mate:')
            print(bamfile.mate(read_dict[i][0]))
    print('=== csn_pair_dict remaining ===')
    if bool(csn_pair_dict):
        for i in csn_pair_dict:
            print(i)
            print(csn_pair_dict[i])
            # print('read remaining:')
            # print(csn_pair_dict[i][0])
            # print('mate:')
            # print(bamfile.mate(csn_pair_dict[i][0]))

    # ===== write tag family size dictionary to file =====
    import pickle
    tag_file = open(args.outfile.split('.sscs')[0] + '.read_families.txt', 'ab+')
    pickle.dump(tag_dict, tag_file)
    tag_file.close()

    summary_stats = '''Total reads overlapping bedfile: {} \n
Unmapped reads: {} \n
Reads with unmapped mate: {} \n
Secondary/Supplementary reads: {} \n
SSCS reads: {} \n
Singletons: {} \n
'''.format(counter, unmapped, unmapped_mate, multiple_mapping, SSCS_reads, singletons)

    stats.write(summary_stats)
    print(summary_stats)

    # === QC ===
    print('===QC metric===')
    print('Mapped reads overlapping bed file should be equivalent with mapped reads in bam file.')
    print('Total mapped reads in bam file: {}'.format(bamfile.mapped))
    # print('Total unmapped reads: {}'.format(bamfile.unmapped))
    # print('Total reads: {}'.format(bamfile.mapped + bamfile.unmapped))

    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    badRead_bam.close()
    tag_bam.close()


    # ===== Create tag family size plot =====
    # Count number of families containing each number of read (e.g. Counter({1: 3737, 32: 660... ->
    # 3737 families are singletons)

    fam_per_read_group = collections.Counter([i for i in tag_dict.values()])
    lst_fam_per_read = list(fam_per_read_group.items()) # conver to list

    total_reads = sum(tag_dict.values())
    # Multiply number of families by read num to get total number of reads in that read group, divide it by total reads
    # to obtain read fraction
    read_fraction = [(i*j)/total_reads for i,j in lst_fam_per_read]

    plt.bar(list(fam_per_read_group), read_fraction)

    # == Determine read family size range to standardize plot axis ==
    plt.xlim([0, math.ceil(lst_fam_per_read[-1][0]/10) * 10])

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

