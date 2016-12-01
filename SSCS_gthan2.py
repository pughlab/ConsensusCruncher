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

from consensus_helper_sc import *

###############################
#          Functions          #
###############################


def coor_separator(chr_lst, chr_len, bedfile=None):
    '''(list, list, str) -> list
    Return list of coordinates based on chromosome arm or bed file (if provided).

    - ChrM not divided by arms
    - Chr arm positions are used to separate bam file into more manageable chunks to avoid running out of memory

    Input:
    - chr_lst
    ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    - chr_len
    [16571, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
    135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520,
    48129895, 51304566, 155270560, 59373566]
    '''
    chr_arm_coor = collections.OrderedDict()

    if bedfile == None:
        # if no bed file provided, use cytoband file instead and separate by chr arm
        filepath = os.path.abspath(inspect.getfile(inspect.currentframe())).rsplit('/', 1)[0]
        bedfile = filepath + '/cytoBand.txt'

    if 'chrM' in chr_lst:
        chr_arm_coor['chrM'] = (0, chr_len[chr_lst.index('chrM')])

    with open(bedfile) as f:
        next(f)  # Skip header
        for line in f:
            chr_arm = line.split('\t')
            chr_key = '{}_{}'.format(chr_arm[0], chr_arm[3])
            start = int(chr_arm[1])
            end = int(chr_arm[2])
            chr_val = (start, end)

            chr_arm_coor[chr_key] = chr_val

            # === HPV Genome === -> Should write script to incorporate any chr not found in cytoband file....
            # if 'chrHPV16_gi_333031' in chr_lst:
            # chr_arm_coor['chrHPV16_gi_333031'] = (0, chr_len[chr_lst.index('chrHPV16_gi_333031')])

    return chr_arm_coor


def insert_pos(cigar):
    '''(list) -> list
    Return 0-based index of positions with insertions, based on list of operations in cigar string.
    [(position, number of inserted bases)]

    === CIGAR ===
    - soft clips, hard clips, insertions, deletions
    - Soft clips are shown in seq, hard clips are not.

    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip

    Examples:
    >>> insert_pos([(4, 14), (0, 29), (2, 2), (0, 30)])
    None
    >>> insert_pos([(0, 47), (2, 3), (0, 27), (1, 6), (0, 13), (4, 30)])
    [(74, 6)]
    >>> insert_pos([(0, 35), (1, 4), (0, 16), (1, 5), (0, 32), (5, 31)])
    [(35, 4), (55, 5)]
    '''
    inserts = []
    index = 0
    for operation in cigar:
        if operation[0] == 1:  # Insertion
            inserts.append((index, operation[1]))
            index += operation[1]
        elif operation[0] == 0:  # Match/mismatch
            index += operation[1]
        else:
            pass

    if not inserts:
        return None

    return sorted(inserts)


def mismatch_pos(cigar, mismatch_tag, readLength):
    '''(list, str) -> lst
    Return 0-based index of mismatch positions in sequence (including insertions, ignoring deletions)
    - list of non-reference allele positions. Positions are adjusted to reflect insertions.

    === Mismatch_tag/MD ===
    - position of mismatch
    ** NOTE: MD only represent information about read aligned to reference, don't take into consideration soft clips,
    insertions, or deletions

    e.g. mismatch_pos([(0, 28), (1, 1), (0, 69)], '19T74G2')
    total bases: cigar = 98, MD = 97

    ========================
    ** Mismatch positions need to account for both cigar and MD (offset MD positions when insertions are present)

    === Test cases ===
    # TTTA_24_58979470_24_58979482_28M5I110M4S_16M5I110M16S_neg_fwd_R2
    # total length 147, soft clip = 3, thus length 143
    # 28M5I110M4S, cigar = 147, MD = 111
    >>> mismatch_pos([(0, 28), (1, 5), (0, 110), (4, 4)], '0T9T23C4C3A9T2T15A3A5C0T0G4T9A0A6A30', )


    # AGTT_24_13843933_24_13843945_8S43M1D5M1I90M_5S115M27S_neg_rev_R1
    # 5S115M27S, cigar = 147, MD = 115
    >>> mismatch_pos([(4, 5), (0, 115), (4, 27)], '20G11A19A10C10T1G1A2G2A1G8G19')


    # AGTT_24_13843945_24_13843933_8S43M1D5M1I90M_5S115M27S_neg_fwd_R2
    # 8S43M1D5M1I90M, cigar = 147, MD = 138 (excluding soft clips, insertions, and deletions), tlen = 139
    >>> mismatch_pos([(4, 8), (0, 43), (2, 1), (0, 5), (1, 1), (0, 90)], '13G29^G0G6C60T15G3A6')






    # 14S29M2D30M50S, cigar = 59, MD = 59
    >>> mismatch_pos([(4, 14), (0, 29), (2, 2), (0, 30), (4, 50)], '13c0g14^ga3g23N0N0N0', 123)
    [14, 27, 28, 46, 70, 71, 72, 73]

    # 47M3D27M6I13M30S
    >>>mismatch_pos([(0, 47), (2, 3), (0, 27), (1, 6), (0, 13), (4, 30)], '19a27^ggc12g25N0N0', 123)
    [0, 19, 59, 74, 75, 76, 77, 78, 79, 91, 92, 93]

    # 35M4I16M5I32M31H <- realistically, hard clips are filtered out
    >>>mismatch_pos([(0, 35), (1, 4), (0, 16), (1, 5), (0, 32), (5, 31)], '3g10g10g1g15c38N0', 92)
    [0, 3, 14, 25, 27, 35, 36, 37, 38, 47, 55, 56, 57, 58, 59, 91, 92]

    8S115M
    >>>mismatch_pos([(4, 8), (0, 115)], '0N0N113', 123)
    [8, 8, 9, 123]

    E.g. mismatch_pos([(0, 19), (2, 1), (0, 79)], '19^A8G70', 98)
    [27]
    ^ = deletion
    19 matches, followed by 1bp deletion, 8 matches, with a G on the reference which is different form aligned read base, and 70 matches

    First mismatch "G' at position 19 + 8
    '''
    # Non-ref allele positions
    mismatches = re.split('[^0-9, \^]+', mismatch_tag)  # split by letters
    mis_pos = []
    index = 0

    inserts = insert_pos(cigar)
    old_inserts = []

    for base in range(len(mismatches)-1):
        if mismatches[base] == '':  # mismatch in first pos
            mis_pos.append(0)
        elif '^' in mismatches[base]:  # deletions are ignored
            index += int(mismatches[base].replace("^", ""))
        else:
            index += int(mismatches[base])
            mis_pos.append(index)
            index += 1  # why add 1??

            # incorporate insertion into index
            # What if insertion in first position???
            for i in inserts:
                if index > i[0] and i not in old_inserts:
                    mis_pos[-1] = mis_pos[-1] + i[1]
                    for val in range(i[1]):
                        mis_pos.append(i[0] + val)
                        index += 1
                    old_inserts.append(i)

    return sorted(mis_pos)


def consensus_maker(readList, readLength, cutoff):
    '''(list, int, int) -> str
    Return consensus sequence and quality score.

    Majority rules concept where if no majority is reached above the cutoff, an 'N' is assigned to the position.
    - Do we want to assign something instead of N?

    - At each position, add quality score to list corresponding to nucleotide.
      Take max quality score of nucleotide with highest frequency

    Bases below the phred quality cutoff (0.7) are excluded from consenus making
    '''
    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    consensus_read = ''
    quality_consensus = []
    proportion_scores = []

    mismatch_pos_lst = []

    for read in readList:
        mismatch_pos_lst.append(mismatch_pos(read.cigartuples, read.get_tag('MD'), read.infer_query_length()))

    for i in range(readLength):
        position_score = [0, 0, 0, 0, 0] # A, C, G, T, N
        quality_score = [[], [], [], [], []]
        phred_fail = 0

        for j in range(len(readList)):
            # === Add N and set qual score to 0 for soft clipped regions ===
            #if i < mismatch_pos_lst[j][0] or i >= mismatch_pos_lst[j][-1]:
                #position_score[4] += 1
                #quality_score[4].append(0)
            #else:
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
            # If there's more than one max, set to Ns
            max_nuc = max_nuc_index[randint(0, len(max_nuc_index)-1)]

            ### In extreme dilutions, Look at real mutations to see if proportion cutoff is setting those real variants into Ns

            ## Make workflow of logic of code

            # === Molecular phred quality of error (multiple probabilities of error) ===
            # error_qualities = quality_score[:max_nuc] + quality_score[(max_nuc +1):]
            ## error qual of majority instead of minority
            err_qual_unlisted = quality_score[max_nuc]
            from itertools import chain
            import math
            #err_qual_unlisted = list(chain(*error_qualities))

            if err_qual_unlisted == []:
                # No error/all consensus, calculate probability of error corresponding to max capped quality score and multiply all of them
                mol_qual = 62
            else:
                P = 1
                for Q in err_qual_unlisted:
                    P *= 10**(-(Q/10))

                mol_qual = round(-10 * math.log10(P))

                if mol_qual > 62:
                    mol_qual = 62

            # frequency of nuc at position > cutoff
            prop_score = position_score[max_nuc]/(len(readList) - phred_fail)
            if prop_score >= cutoff:
                consensus_read += nuc_lst[max_nuc]
                quality_consensus.append(mol_qual)
                proportion_scores.append(prop_score)
            else:
                raise ValueError

        except:
            consensus_read += 'N'
            quality_consensus.append(0)
            proportion_scores.append(0)

    return consensus_read, quality_consensus, proportion_scores
# , position_score[max_nuc], phred_fail # INCLUDE PHRED FAIL AS THIS DIFFERS FROM TOTAL NUMBER OF READS IN FAMILY


def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--cutoff", action="store", dest="cutoff", help="nucleotide base % cutoff", required=True)
    # parser.add_argument("--Ncutoff", action = "store", dest="Ncutoff", help="N % cutoff", required = True)
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output SSCS BAM file", required=True)
    parser.add_argument("--bedfile", action="store", dest="bedfile", help="input bedfile for data division", required=False)
    parser.add_argument("--taginRGbam", action="store", dest="RGbam",
                        help="output bamfile with consensus tag in read group", required=False)
    args = parser.parse_args()

    start_time = time.time()

    # ===== Initialize input and output bam files =====
    bamfile = pysam.AlignmentFile(args.infile, "rb")
    SSCS_bam = pysam.AlignmentFile(args.outfile, "wb", template = bamfile)
    stats = open('{}.stats.txt'.format(args.outfile.split('.sscs')[0]), 'w')
    singleton_bam = pysam.AlignmentFile('{}.singleton.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)
    badRead_bam = pysam.AlignmentFile('{}.badReads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    tag_bam = pysam.AlignmentFile('{}.tag_reads.bam'.format(args.outfile.split('.sscs')[0]), "wb", template = bamfile)

    # set up fastq files
    fastqFile1 = open('{}.sscs_R1.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')
    fastqFile2 = open('{}.sscs_R2.fastq.gz'.format(args.outfile.split('.sscs')[0]), 'w')

    # set up time tracker
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.sscs')[0]), 'w')

    # ===== Initialize dictionaries =====
    read_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    pair_dict = collections.defaultdict(list)

    quality_dict = collections.defaultdict(list)
    prop_dict = collections.defaultdict(list)

    # ===== Initialize counters =====
    unmapped = 0
    unmapped_flag = 0
    bad_reads = 0  # secondary/supplementary reads
    poor_mapq = 0
    counter = 0
    singletons = 0
    SSCS_reads = 0

    # ===== Determine data division coordinates =====
    # division by bed file if provided, or by chr arm coordinates
    if 'args.bedfile' in locals():
        division_coor = coor_separator(chrm, chr_len, args.bedfile)
    else:
        chrm = [x['SN'] for x in bamfile.header['SQ']]
        chr_len = [x['LN'] for x in bamfile.header['SQ']]
        division_coor = coor_separator(chrm, chr_len)

    # ===== Process data in chunks =====
    for x in division_coor.keys():
        chr_data = read_bam(bamfile,
                            pair_dict = pair_dict,
                            read_dict = read_dict,
                            tag_dict = tag_dict,
                            badRead_bam = badRead_bam,
                            read_chr = x.rsplit('_', 1)[0],
                            read_start = chr_arm_coor[x][0],
                            read_end = chr_arm_coor[x][1]
                            )

        read_dict = chr_data[0]
        tag_dict = chr_data[1]
        pair_dict = chr_data[2]

        counter += chr_data[3]
        unmapped += chr_data[4]
        unmapped_flag += chr_data[5]
        bad_reads += chr_data[6]
        poor_mapq += chr_data[7]

        # ===== Determine read length =====
        # Randomly picked family with largest family size (most PCR dupes)
        max_family_size = max(tag_dict.values())
        family = [key for key in tag_dict.items() if key[1] == max_family_size]
        rand_read = read_dict[choice(list(family))[0]]

        # Infer length from most common cigar string
        cigar_mode = statistics.mode([r.cigarstring for r in rand_read])
        cigar_mode_read = [r for r in rand_read if r.cigarstring == cigar_mode][0]
        readLength = cigar_mode_read.infer_query_length()

        # ===== Create consensus sequences as paired reads =====
        written_pairs = []

        for readPair in list(pair_dict.keys()):
            if len(pair_dict[readPair]) == 2:
                for tag in pair_dict[readPair]:
                    # === Check for singletons ===
                    if tag_dict[tag] < 2:
                        singletons += 1
                        singleton_bam.write(read_dict[tag][0])
                    else:
                        SSCS = consensus_maker(read_dict[tag], readLength, float(args.cutoff))

                        SSCS_read = create_aligned_segment(read_dict[tag], SSCS[0], SSCS[1])

                        query_name = readPair + ':' + str(tag_dict[tag])
                        SSCS_read.query_name = query_name

                        #quality_dict[query_name] += [SSCS[1]]
                        #prop_dict[query_name] += [SSCS[2]] #### SAVE BY QUERY NAME OR TAG NAME?
                        #tag_quality_dict[tag_dict[tag]] += [round(np.mean(SSCS[1]))]

                        #try:
                        #SSCS_read.set_tag('PR', SSCS[2])
                        #except:
                            #print(tag)
                            #print(SSCS)
                            #print(SSCS[2])

                        #aligned_seq = SSCS_read.query_alignment_sequence
                        #if aligned_seq.count('N')/len(aligned_seq) > float(args.Ncutoff):
                            #print('uh oh too many Ns')
                            #print(aligned_seq)
                            #print(SSCS_read)
                            #continue

                        # === Write new bamfile with consensus query name in read group ===
                        if 'args.RGbam' in locals():
                            for r in read_dict[tag]:
                                r.set_tag('RG', query_name)
                                tag_bam.write(r)

                        # ===== Write consensus bam =====
                        SSCS_bam.write(SSCS_read)
                        SSCS_reads += 1

                        # ===== write as fastq file =====
                        #if SSCS_read.is_reverse and SSCS_read.is_read1:
                        #fastq_seq = SSCS_read.query_sequence.decode("utf-8")
                        fastq_seq = SSCS[0]

                        if 'rev' in tag:
                            fastq_seq = reverse_seq(fastq_seq)
                            fastq_qual = pysam.qualities_to_qualitystring(reversed(SSCS[1]))
                        else:
                            fastq_qual = pysam.qualities_to_qualitystring(SSCS[1])

                        if 'R1' in tag:
                            fastqFile1.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))
                        else:
                            fastqFile2.write('@{}\n{}\n+\n{}\n'.format(SSCS_read.query_name, fastq_seq, fastq_qual))

                    # Remove read from dictionary after writing
                    del read_dict[tag]

                # Remove key from dictionary after writing
                del pair_dict[readPair]

                ##### NEED TO ALSO REMOVE SC FROM SC_LST AFTER WRITING PAIR!!!!
        print(x)

        # This ends up writing everything to the file, also dictionary of reads is not decreasing over time -> memory issue
        #read_dict_file = open(args.outfile.split('.sscs')[0] + '.read_dict.txt', 'ab+')
        #pickle.dump(bam_dict, read_dict_file)
        #read_dict_file.close()


        print(str((time.time() - start_time)/60))
        try:
            time_tracker.write(x + ': ')
            time_tracker.write(str((time.time() - start_time)/60) + '\n')

        except:
            continue


    # Reads left over at the end due to how pysam fetches reads -> does it get read twice?
    print('read pairs remaining')
    if bool(read_pair_dict):
        for i in read_pair_dict:
            print(i)
            print(read_pair_dict[i])
            print(bamfile.mate(read_pair_dict[i][0]))

    print('pair_dict remaining')
    if bool(pair_dict):
        for i in pair_dict:
            print(i)
            print(pair_dict[i])
            print(read_dict[pair_dict[i][0]])
            print(bamfile.mate(read_dict[pair_dict[i][0]][0]))
    print('bam_dict remaining')
    if bool(read_dict):
        for i in read_dict:
            print(i)
            print(read_dict[i])
            print(read_dict[i][0])
            #print(bamfile.mate(bam_dict[i][0]))

    # ===== write tag family size dictionary to file =====
    import pickle
    #qual_file = open(args.outfile.split('.sscs')[0] + '.q_scores.txt', 'ab+')
    #pickle.dump(quality_dict, qual_file)
    #qual_file.close()
    #quality_dict = collections.defaultdict(list)

    #prop_file = open(args.outfile.split('.sscs')[0] + '.prop_scores.txt', 'ab+')
    #pickle.dump(prop_dict, prop_file)
    #prop_file.close()
    #prop_dict = collections.defaultdict(list)

    # (key = tags, value = int [number of reads in that family])
    tag_file = open(args.outfile.split('.sscs')[0] + '.read_families.txt', 'ab+')
    pickle.dump(tag_dict, tag_file)
    tag_file.close()


    summary_stats = '''Total reads: {} \n
Unmapped reads: {} \n
Unmapped flag reads: {} \n
Secondary/Supplementary reads: {} \n
Poor mapping quaity (<5): {} \n
SSCS reads: {} \n
Singletons: {} \n
'''.format(counter, unmapped, unmapped_flag, bad_reads, poor_mapq, SSCS_reads, singletons)

    stats.write(summary_stats)

    time_tracker.close()
    stats.close()
    bamfile.close()
    SSCS_bam.close()
    badRead_bam.close()
    tag_bam.close()
    fastqFile1.close()
    fastqFile2.close()

    # ===== Create tag family size plot =====
    # Count number of families containing each number of read (e.g. Counter({1: 3737, 32: 660... -> 3737 families are singletons)

    fam_per_read_group = collections.Counter([i for i in tag_dict.values()])
    lst_fam_per_read = list(fam_per_read_group.items()) # conver to list

    total_reads = sum(tag_dict.values())
    # Multiply number of families by read num to get total number of reads in that read group, divide it by total reads to obtain read fraction
    read_fraction = [(i*j)/total_reads for i,j in lst_fam_per_read]

    plt.bar(list(fam_per_read_group), read_fraction)

    # == Determine read family size range to standardize plot axis ==
    plt.xlim([0, math.ceil(lst_fam_per_read[-1][0]/10) * 10])

    #plt.locator_params(axis = 'x', nbins=lst_fam_per_read[-1][0]//5)
    plt.xlabel('Tag family size (# of reads per family)')
    plt.ylabel('Fraction of total reads')

    plt.savefig(args.outfile.split('.sscs')[0]+'_tag_fam_size.png')

    # ===== Create prop plot =====



###############################
##           Main            ##
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    print((time.time() - start_time)/60)

