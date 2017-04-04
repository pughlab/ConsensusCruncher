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
import matplotlib.pyplot as plt
import math
import time

from consensus_helper import *

###############################
#       Helper Functions      #
###############################


def genomicBasedCigar(cigar, pos):
    '''(str, int) -> list
    Return list of genomic positions based on cigar and start coordinate.

    Function accounts for insertions in genome coordinate.
    '''
    pattern = re.compile('([MIDNS=])')
    C = pattern.split(cigar)[1::2]  # Cigar e.g. ['S', 'M']
    I = pattern.split(cigar)[::2][:-1]  # Index e.g. ['33', '90']
    cord = []
    for i in range(len(C)):
        if C[i] == 'S':
            cord = cord + [0] * int(I[i])
        elif C[i] == 'M':
            cord = cord + list(range(pos, pos + int(I[i]), 1))
            pos += int(I[i])
        elif C[i] == 'I':
            cord = cord + [pos] * int(I[i])
        elif C[i] == 'D':
            pos += int(I[i])

    return cord


def consensus_maker(readList, cutoff, failed_bases):
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
    # === Genome coordinate of each base ===
    # As each family shares the same start coordinate and cigar, this info can be collected from any read
    base_coor = genomicBasedCigar(readList[0].cigarstring, readList[0].reference_start)

    # === Determine read length ===
    cigar_mode = collections.Counter([r.cigarstring for r in readList]).most_common(1)[0][0]
    cigar_mode_read = [r for r in readList if r.cigarstring == cigar_mode][0]
    readLength = cigar_mode_read.infer_query_length()

    # === Initialize counters ===
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = []
    proportion_scores = []

    for i in range(readLength):
        nuc_count = [0, 0, 0, 0, 0]  # A, C, G, T, N
        failed_nuc_count = [0, 0, 0, 0, 0]
        quality_score = [[], [], [], [], []]
        phred_fail = 0

        for j in range(len(readList)):
            # === Phred filter Q30 cut-off ===
            if readList[j].query_qualities[i] < 30:
                nuc = readList[j].query_sequence[i]
                nuc_index = nuc_lst.index(nuc)
                failed_nuc_count[nuc_index] += 1
                # quality_score[4] += 0
                phred_fail += 1
            else:
                nuc = readList[j].query_sequence[i]
                nuc_index = nuc_lst.index(nuc)
                nuc_count[nuc_index] += 1
                quality_score[nuc_index].append(readList[j].query_qualities[i])

        # Find most frequent nuc (don't worry about ties (2 maxes) as it won't pass proportion cut-off and N will
        # be assigned)
        max_nuc_index = nuc_count.index(max(nuc_count))

        # === Molecular phred quality (consensus quality score) ===
        max_nuc_quality = quality_score[max_nuc_index]

        base_fail = False

        if max_nuc_quality != []:
            P = 1
            for Q in max_nuc_quality:
                P *= 10**(-(Q/10))

            # large families leads to multiplication of numerous small floats resulting in zero
            if P == 0:
                mol_qual = 62
            else:
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
            prop_score = nuc_count[max_nuc_index]/phred_pass_reads
            if prop_score >= cutoff:
                if base_fail == True:
                    # test to see if a position that has failed quality score making can still pass filters
                    print('base fail == True!!')
                consensus_read += nuc_lst[max_nuc_index]
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

            # === Write failed bases to file ===
            # position of 2 most freq bases
            if nuc_count == [0, 0, 0, 0, 0]:
                max_nuc_index = failed_nuc_count.index(max(failed_nuc_count))
                max_nuc = '{}*'.format(nuc_lst[max_nuc_index])
                failed_nuc_count[max_nuc_index] = 0
                if failed_nuc_count == [0, 0, 0, 0, 0]:
                    second_max_nuc = 'NA'
                else:
                    second_max_index = failed_nuc_count.index(max(failed_nuc_count))
                    second_max_nuc = '{}*'.format(nuc_lst[second_max_index])
            else:
                max_nuc = nuc_lst[max_nuc_index]
                nuc_count[max_nuc_index] = 0
                if nuc_count == [0, 0, 0, 0, 0]:
                    second_max_nuc = 'NA'
                else:
                    second_max_index = nuc_count.index(max(nuc_count))
                    second_max_nuc = nuc_lst[second_max_index]

            # chr number
            if readList[0].reference_id == 0:
                chr = 'M'
            elif readList[0].reference_id == 23:
                chr = 'X'
            elif readList[0].reference_id == 24:
                chr = 'Y'
            else:
                chr = readList[0].reference_id
            # print(nuc_count)
            # print('max nuc {}'.format(max_nuc))
            # print('2nd max nuc {}'.format(second_max_nuc))
            failed_info = 'chr{}\t{}\t{}\t{}\n'.format(chr,
                                                       base_coor[i],  # Pos
                                                       max_nuc,  # Most freq base
                                                       second_max_nuc)  # 2nd most freq base

            failed_bases.write(failed_info)

    return consensus_read, quality_consensus, proportion_scores


###############################
#        Main Function        #
###############################

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
    SSCS_uncollapesd = pysam.AlignmentFile('{}.sscs.uncollapsed.bam'.format(args.outfile.split('.sscs')[0]), 'wb', template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    failed_bases = open('{}.failed_bases.txt'.format(args.outfile.split('.sscs')[0]), 'w')

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
    SSCS_uncollapsed_reads = 0

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
                        # === Write uncollapsed SSCS reads ===
                        # write reads from tag family sizes >= 2 to file prior to collapsing
                        for read in read_dict[tag]:
                            SSCS_uncollapesd.write(read)
                            SSCS_uncollapsed_reads += 1

                        # === Create collapsed SSCSs ===
                        SSCS = consensus_maker(read_dict[tag], float(args.cutoff), failed_bases)

                        query_name = readPair + ':' + str(tag_dict[tag])
                        SSCS_read = create_aligned_segment(read_dict[tag], SSCS[0], SSCS[1], query_name)

                        # === Write consensus bam ===
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
            try:
                print(i)
                print('read remaining:')
                print(pair_dict[i][0])
                print('mate:')
                print(bamfile.mate(pair_dict[i][0]))
            except ValueError:
                print("Mate not found")
    print('=== read_dict remaining ===')
    if bool(read_dict):
        for i in read_dict:
            try:
                print(i)
                print('read remaining:')
                print(read_dict[i][0])
                print('mate:')
                print(bamfile.mate(read_dict[i][0]))
            except ValueError:
                print("Mate not found")
    print('=== csn_pair_dict remaining ===')
    if bool(csn_pair_dict):
        for i in csn_pair_dict:
            try:
                print(i)
                print(csn_pair_dict[i])
                # print('read remaining:')
                # print(csn_pair_dict[i][0])
                # print('mate:')
                # print(bamfile.mate(csn_pair_dict[i][0]))
            except ValueError:
                print("Mate not found")


    # ===== write tag family size dictionary to file =====
    # count of tags within each family size
    tags_per_fam_size = collections.Counter([i for i in tag_dict.values()])
    lst_fam_per_read = list(tags_per_fam_size.items())  # convert to list
    with open(args.outfile.split('.sscs')[0] + '.read_families.txt', "w") as stat_file:
        stat_file.write('family_size\tfrequency\n')
        stat_file.write('\n'.join('%s\t%s' % x for x in lst_fam_per_read))

	# == Pandas method ==
#     import pandas as pd
#     tag_df = pd.DataFrame(list(tag_dict.items()), columns=['tag_ID', 'family_size'])
#     tag_df_summary = tag_df.join(tag_df['tag_ID'].str.split('_', expand=True))
#     tag_df_summary.columns = ['tag_ID', 'family_size', 'barcode', 'R1chr', 'R1start', 'R2chr', 'R2start', 'R1cigar',
#                               'R2cigar', 'strand', 'orientation', 'RG']
#     tag_df_summary.to_csv(args.outfile.split('.sscs')[0] + '.read_families.txt', index=None, sep='\t', mode='a')

	# == Pickle method ==
    # import pickle
    # tag_file = open(args.outfile.split('.sscs')[0] + '.read_families.txt', 'ab+')
    # pickle.dump(tag_dict, tag_file)
    # tag_file.close()

    # === STATS ===
    # Note: total reads = unmapped + secondary + SSCS uncollapsed + singletons
    summary_stats = '''Original - Total reads overlapping bedfile: {}
Original - Unmapped reads: {}
Original - Secondary/Supplementary reads: {}
SSCS reads: {}
SSCS uncollapsed reads: {}
Singletons: {} \n'''.format(counter, unmapped, multiple_mapping, SSCS_reads, SSCS_uncollapsed_reads, singletons)

    stats.write(summary_stats)
    print(summary_stats)

    # === QC ===
    # print('===QC metric===')
    print('# QC: Mapped reads overlapping bed file should be equivalent to mapped reads in bam file.')
    print('Total mapped reads in bam file: {}'.format(bamfile.mapped))
    # print('Total unmapped reads: {}'.format(bamfile.unmapped))
    # print('Total reads: {}'.format(bamfile.mapped + bamfile.unmapped))

    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    SSCS_uncollapesd.close()
    badRead_bam.close()
    failed_bases.close()

    # ===== Create tag family size plot =====
    # Count number of families containing each number of read (e.g. Counter({1: 3737, 32: 660... ->
    # 3737 families are singletons)

    fam_per_read_group = collections.Counter([i for i in tag_dict.values()])
    lst_fam_per_read = list(fam_per_read_group.items())  # convert to list

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
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)

