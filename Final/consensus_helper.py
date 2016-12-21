#!/usr/bin/env python

###############################################################
#
#                          Consensus Helper
#
# Author: Nina Wang
# Date Created: Mar 24, 2016
###############################################################
# Function:
# Written for Python 3.5.1
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
import os
import inspect


###############################
#          Functions          #
###############################


def bed_separator(bedfile):
    '''(str) -> dict
    Return list of coordinates based on bedfile.
    '''
    coor = collections.OrderedDict()

    with open(bedfile) as f:
        for line in f:
            chr_arm = line.split('\t')
            chr_key = '{}_{}'.format(chr_arm[0], chr_arm[3])
            start = int(chr_arm[1])
            end = int(chr_arm[2])
            chr_val = (start, end)

            coor[chr_key] = chr_val

    return coor


def which_read(flag):
    '''(int) -> str
    Returns read 1 or 2 based on flag.

    Test cases:
    >>> which_read(83)
    'R1'
    >>> which_read(131)
    'R2'
    >>> which_read(177)
    'R2'
    '''
    read1 = [99, 83, 67, 115, 81, 97, 65, 113]
    read2 = [147, 163, 131, 179, 161, 145, 129, 177]

    if flag in read1:
        read = 'R1'
    elif flag in read2:
        read = 'R2'
    else:
        print('UNMAPPED READ ERROR')
        print(flag)
        read = None

    return read


def which_strand(read):
    '''(pysam.calignedsegment.AlignedSegment) -> str
    Return DNA strand based on flags.

    Notes: Need to determine orientation for translocations (flags 65, 129, 113, 177) where pair occurs in same direction
    - pos: read is first in pair (R1) and has lower chr number than mate with higher chr number (read_chr < mate_chr) OR
           read is second in pair (R2) and has higher chr number than mate with lower chr number (read_chr > mate_chr)
           e.g. H1080:278:C8RE3ACXX:6:1209:19123:36119|CTCT	113	12	25398309	60	98M	16	7320258	98 -> (chr12_chr16_R1) 'pos'
                H1080:278:C8RE3ACXX:6:1209:19123:36119|CTCT	177	16	7320258	60	98M	12	25398309	98 -> (chr16_chr12_R2) 'pos'
    - neg: Opposite cases of pos
           e.g. H1080:278:C8RE3ACXX:6:1307:20616:55254|CTCT	177	12	25398309	60	98M	16	7320258	98 -> (chr16_chr12_R1) 'neg'
                H1080:278:C8RE3ACXX:6:1307:20616:55254|CTCT	113	16	7320258	60	98M	12	25398309	98 -> (chr12_chr16_R2) 'neg'

    Test cases:
    Flag = 147 -> 'pos'
    Flag = 67 -> 'pos'
    Flag = 131 -> 'pos'
    Flag = 81 -> 'neg'
    '''
    # Flags indicating strand direction
    pos = [99, 147, 67, 131, 97, 145]  # note 65/129 direction not defined
    neg = [83, 163, 115, 179, 81, 161]  # note 113/177 direction not defined
    no_ori = [65, 129, 113, 177]

    if read.flag in pos:
        strand = 'pos'
    elif read.flag in neg:
        strand = 'neg'
    elif read.flag in no_ori:
        if (read.reference_id < read.next_reference_id and which_read(read.flag) == 'R1') or \
                (read.reference_id > read.next_reference_id and which_read(read.flag) == 'R2'):
            strand = 'pos'
        else:
            strand = 'neg'
    else:
        # Only uniquely mapped reads (with flags indicated above) should be retained as 'bad reads' were filtered out at
        # a previous step
        print('STRAND ERROR')
        print(read.flag)
        strand = None

    return strand


def cigar_order(read, read_pair):
    '''(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> str
    Return ordered cigar string from R1 and R2

    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip

    Pos strand, R1 cigar first
    Neg strand, R2 cigar first
    => want cigar in same order for Duplex consensus formation

    Examples:
        (+) [99]  137M10S    => 137M10S_147M
            [147] 147M

        (-) [83]  147M   => 137M10S_147M
            [163] 137M10S
    '''
    ori_strand = which_strand(read)
    read_num = which_read(read.flag)

    if (ori_strand == 'pos' and read_num == 'R1') or \
            (ori_strand == 'neg' and read_num == 'R2'):
        cigar = '{}_{}'.format(read.cigarstring,
                               read_pair.cigarstring)
    else:
        cigar = '{}_{}'.format(read_pair.cigarstring,
                               read.cigarstring)

    return cigar


def unique_tag(read, cigar, barcode):
    '''(pysam.calignedsegment.AlignedSegment, str, str) -> str
    Return unique identifier tag for one read of a strand of an individual molecule.

    Tag uses following characteristics to group reads belonging to the same strand of an individual molecule (PCR dupes):

    barcode_ReadChr_ReadStart_MateChr_MateStart_cigar_strand_orientation_readNum
    TTTG_24_58847416_24_58847448_137M10S_147M_pos_fwd_R1

    Notes:
        - paired reads have cigar in same order (see cigar_order fx for details)
        - barcode location differs depending on whether you're making SSCS vs DCS
        - chr of R1 and R2 included as they might differ in case of translocation
        - orientation included for cases where mate and complimentary strand are the same (e.g. flags 65/129 -> +/+/+/+)

    Example:
       R1 --->   <--- R2
    (+) TT-----------GT

    (-) TT-----------GT
       R2 --->   <--- R1

    ### Change to duplex read example!!! ###

    R1 of (+) -> TTTG_24_58847416_24_58847448_137M10S_147M_pos_fwd_R1
    NB500964:12:HTTG2BGXX:4:22601:26270:1144|TTTG	99	24	58847416	17	137M10S	24	58847448	137

    R2 of (+) -> TTTG_24_58847448_24_58847416_137M10S_147M_pos_rev_R2
    NB500964:12:HTTG2BGXX:4:22601:26270:1144|TTTG	147	24	58847448	17	147M	24	58847416	147

    R1 of (-) -> TGTT_24_58847448_24_58847416_137M10S_147M_neg_rev_R1

    R2 of (-) -> TGTT_24_58847416_24_58847448_137M10S_147M_neg_fwd_R2
    '''

    strand = which_strand(read)

    orientation = 'fwd'
    if read.is_reverse:
        orientation = 'rev'

    readNum = which_read(read.flag)

    # Add flag direction

    # Unique identifier for strand of individual molecules
    tag = '{}_{}_{}_{}_{}_{}_{}_{}_{}'.format(barcode,  # mol barcode
                                              read.reference_id,  # chr num
                                              read.reference_start,  # start R1 (0-based)
                                              read.next_reference_id,
                                              read.next_reference_start,  # start R2
                                              cigar,
                                              strand,
                                              orientation,  # strand direction
                                              readNum
                                              )
    return tag


def sscs_qname(tag, flag):
    '''(str, int) -> str
    Return new query name for consensus sequence: barcode_chr_start_chr_end_strand_flags

    * Since multiple reads go into making a consensus, a new query name is needed as an identifier for consensus read
    pairs. *
    - Read pairs share the same query name to indicate they're mates

    original query name (H1080:278:C8RE3ACXX:6:1308:18882:18072) ->
    new query name (TTTG_24_58847448_24_58847416_137M10S_147M_pos_99_147)

    Reason why we want both the flag and strand:
    - flag: to differentiate between reads with pair in same orientation (e.g. 65 / 129)
    - strand: needed for making duplex consensus sequences as reads are grouped based on strand and not specific flags

    Note: coordinates are ordered from low -> high, so read pairs share the same coordinate

    ISSUE -> pysam converts coordinates to be 0-based, query name coordinates don't match coordinate seen in read

    Examples:
    (+)                                                 [Flag]
    TTTG_24_58847416_24_58847448_137M10S_147M_pos_fwd_R1 [99] --> TTTG_24_58847416_24_58847448_137M10S_147M_pos_99_147
    TTTG_24_58847448_24_58847416_137M10S_147M_pos_rev_R2 [147]

    (-)
    TGTT_24_58847448_24_58847416_137M10S_147M_neg_rev_R1 [83] --> TGTT_24_58847416_24_58847448_137M10S_147M_neg_83_163
    TGTT_24_58847416_24_58847448_137M10S_147M_neg_fwd_R2 [163]

    Special case (mate and complimentary reads are all in the same direction)
    - Use coordinate and flags to differentiate between strand

    Test cases:
    >>> sscs_qname('TTTG_24_58847416_24_58847448_137M10S_147M_pos_fwd_R1', 99)
    'TTTG_24_58847416_24_58847448_137M10S_147M_pos_99_147'
    >>> sscs_qname('TTTG_24_58847448_24_58847416_137M10S_147M_pos_rev_R2', 147)
    'TTTG_24_58847416_24_58847448_137M10S_147M_pos_99_147'
    >>> sscs_qname('TGTT_24_58847448_24_58847416_137M10S_147M_neg_rev_R1', 83)
    'TGTT_24_58847416_24_58847448_137M10S_147M_neg_83_163'
    >>> sscs_qname('TGTT_24_58847416_24_58847448_137M10S_147M_neg_fwd_R2', 163)
    'TGTT_24_58847416_24_58847448_137M10S_147M_neg_83_163'
    '''
    flag_pairings = {99:147, 147:99, 83:163, 163:83,
                     # mapped within insert size, but wrong orientation (++, --)
                     67:131, 131:67, 115:179, 179:115,
                     # mapped uniquely, but wrong insert size
                     81:161, 161:81, 97:145, 145:97,
                     # wrong insert size and wrong orientation
                     65:129, 129:65, 113:177, 177:113
                     }

    ref_chr = tag.split("_")[1]
    mate_chr = tag.split("_")[3]
    ref_coor = tag.split("_")[2]
    mate_coor = tag.split("_")[4]

    # === Order qname by chr coordinate ===
    # e.g. CCCC_12_25398156_12_25398064 -> CCCC_12_25398064_12_25398156
    if (int(ref_coor) > int(mate_coor) and ref_chr == mate_chr) or \
       (int(ref_chr) > int(mate_chr) and ref_chr != mate_coor):
        new_tag = tag.split("_")
        new_tag[1] = mate_chr # swap ref with mate
        new_tag[3] = ref_chr
        new_tag[2] = mate_coor
        new_tag[4] = ref_coor
        new_tag = "_".join(new_tag)[:-7]

    else:
        # for reads with the same coordinate mate (start and stop are the same)
        new_tag = tag[:-7]

    # === Add flag information to query name ===
    # Flags allow further differentiation as read and its duplex might have same coordinate and strand direction
    # smaller flag ordered first
    #
    # if flag < flag_pairings[flag]:
    #     new_tag = '{}_{}_{}'.format(new_tag, flag, flag_pairings[flag])
    # else:
    #     new_tag = '{}_{}_{}'.format(new_tag, flag_pairings[flag], flag)

    return new_tag


def read_bam(bamfile, pair_dict, read_dict, csn_pair_dict, tag_dict, badRead_bam, read_chr = None,
             read_start = None, read_end = None, duplex = None):
    '''(bamfile object, dict, dict, dict, str, int, int, str) ->
    dict, dict, dict, dict, dict, int, int, int, int

    === Input ===
    - bamfile: pysam.AlignmentFile object of bamfile

    - pair_dict: dictionary of paired reads based on query name to process data in pairs (note: values are removed once
                 pair assigned to corresponding dict)
                 -> retains data from translocations or reads crossing bam division regions (e.g. cytobands) to preserve
                    pairing as file is divided into sections
                {query name: [read 1, read 2]}

    - read_dict: dictionary of reads sharing a common tag (aka from the same unique molecule)
                 {read_tag: [<pysam.calignedsegment.AlignedSegment>, <pysam.calignedsegment.AlignedSegment>, ..etc.]}

    - csn_pair_dict: dictionary of paired tags based on consensus tag
                     -> reads are removed from pair_dict and added to read_dict & csn_pair_dict to track pairs
                     {consensus_tag: [R1_tag, R2_tag]}

    - tag_dict: integer dictionary indicating number of reads in each read family
                 {read_tag: 2, ..etc}

    - badRead_bam: bamfile recording bad reads

    -- Optional --
    # For large bamfiles that are split into regions
    - read_chr: string of chromosome region to fetch reads
    - read_start: integer of starting position to fetch reads
    - read_end: integer of stopping position to fetch reads

    # For duplex consensus makingn
    - duplex: any string or bool specifying duplex consensus making [e.g. TRUE], necessary for parsing data as
              query name for SSCS and DCS differ

    === Output ===
    1) read_dict: dictionary of bamfile reads
                 {read_tag: [<pysam.calignedsegment.AlignedSegment>,
                 <pysam.calignedsegment.AlignedSegment>, ..etc.]}
                 - Key: barcode_chr_startR1_startR2_strand_ReadNum
                 - Value: list of bamfile reads

    2) tag_dict: integer dictionary indicating number of reads in each read family
                 {read_tag: 2, ..etc}

    3) pair_dict: dictionary of paired tags
                    {query name: [read_tag, mate_tag]}

    4) counter: total number of reads

    5) unmapped: reads without a flag

    6) unmapped_flag: unmapped reads or unmapped mates

    7) bad_reads: number of reads that not properly mapped
                  - secondary reads: same sequence aligns to multiple locations
                  - supplementary reads: multiple parts of sequence align to
                                         multiple locations
    8) poor_mapq: number of poorly mapped reads (mapping quality < 5)
    '''
    # ===== Fetch data given coordinates =====
    if read_chr == None:
        bamLines = bamfile.fetch(until_eof = True)
    else:
        bamLines = bamfile.fetch(read_chr, read_start, read_end)

    # ===== Initialize counters =====
    unmapped = 0
    bad_reads = 0  # secondary/supplementary reads
    counter = 0

    for line in bamLines:
        # Parse out reads that don't fall within region
        if read_chr is not None:
            # this excludes certain reads that were previously seen
            if line.reference_start < read_start or line.reference_start > read_end:
                continue

        counter += 1

        # ===== Filter out 'bad' reads =====
        # filter out reads by pairs -> some mapped reads also filtered as their mate might be unmapped
        bad_flags = [73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141]  # Unmapped reads with mapped mate
        # add statement this does not count unmapped pairs
        badRead = True

        if line.flag in bad_flags:
            unmapped += 1
        elif line.is_secondary:
            bad_reads += 1
        elif line.is_supplementary:
            bad_reads += 1
        else:
            badRead = False

        # Write bad reads to file
        if badRead:
            badRead_bam.write(line)
        else:
            pair_dict[line.qname].append(line)

            # ===== Create tag once reads are paired =====
            if len(pair_dict[line.qname]) == 2:
                cigar = cigar_order(pair_dict[line.qname][0], pair_dict[line.qname][1])

                for read in pair_dict[line.qname]:
                    # ===== Extract molecular barcode =====
                    # Barcodes in diff position for duplex consensus formation
                    if duplex == None or duplex == False:
                        # SSCS query name: H1080:278:C8RE3ACXX:6:1308:18882:18072|CACT
                        barcode = read.qname.split("|")[1]
                    else:
                        # DCS query name: CCTG_12_25398000_12_25398118_neg_83_163:5
                        barcode = read.qname.split("_")[0]

                    # Unique identifier for grouping reads belonging to the same strand of an individual molecule
                    tag = unique_tag(read, cigar, barcode)

                    # Create consensus tag to be used as new query name for consensus reads
                    consensus_tag = sscs_qname(tag, int(read.flag))

                    # ===== Add read pair to dictionary =====
                    if tag not in read_dict and tag not in tag_dict:
                        # New read family
                        read_dict[tag] = [read]
                        tag_dict[tag] += 1
                        # === Track read pair by tag instead of query name in pair_dict ===
                        # Reads grouped by query name are removed from pair_dict once they're added to read_dict,
                        # track paired tags using 'consensus' tags (PCR dupes share same tag)
                        # {qname: [read1, read2]} -> {consensus_tag: [R1_tag, R2_tag]}
                        if consensus_tag not in csn_pair_dict:
                            csn_pair_dict[consensus_tag] = [tag]
                        elif len(csn_pair_dict[consensus_tag]) == 2:
                            # in case multiple reads share same consensus tag (e.g. overlapped reads different by 1 base)
                            consensus_tag += '2'
                            if consensus_tag not in csn_pair_dict:
                                csn_pair_dict[consensus_tag] = [tag]
                            else:
                                csn_pair_dict[consensus_tag].append(tag)
                        else:
                            csn_pair_dict[consensus_tag].append(tag)
                    elif tag in tag_dict and read not in read_dict[tag]:
                        # Reads belonging to the same family as another read (PCR dupes)
                        read_dict[tag].append(read)
                        tag_dict[tag] += 1
                    else:
                        # === Data fetch error ===
                        # If tag found in tag_dict, but not in read_dict means line was previously read and written
                        # to file (tag removed from read_dict once written)
                        print('Pair already written: line read twice - double check to see if its overlapping / near cytoband region (point of data division)')
                        # fetched data region - pysam fetch fx will get read twice if it overlaps fetch coordinates
                        print(read_chr, read_start, read_end)
                        print(tag)
                        print(read)
                        print(tag_dict[tag])
                        print(read_dict[tag][0])
                        print(consensus_tag)
                        print(tag in tag_dict)
                        print(csn_pair_dict[consensus_tag])

                        print(pair_dict[line.qname][0])
                        print(pair_dict[line.qname][1])
                        counter -= 1

                # remove read pair qname from pair_dict once reads added to read_dict as pairs
                pair_dict.pop(line.qname)

    return read_dict, tag_dict, pair_dict, csn_pair_dict, counter, unmapped, bad_reads


def read_mode(field, bam_reads):
    '''(str, lst) -> str
    Return mode (most common occurrence) of specified field
    e.g. cigarstring, flag, mapping quality, template_length
    '''
    field = 'i.{}'.format(field)
    # Rank by number of occurrences
    field_lst = collections.Counter(eval(field) for i in bam_reads).most_common()
    # Take max occurrences
    common_field_lst = [i for i, j in field_lst if j == field_lst[0][1]]
    # Randomly select max if there's multiple
    common_field = common_field_lst[randint(0, len(common_field_lst)-1)]

    return common_field


def consensus_flag(bam_reads):
    '''(list) -> str
    Return consensus flag given list of reads from the same family.

    If multiple flags are present within same family and a max can't be determined, prioritize flags in correct orientation
    and within insert size.
    e.g.
    H1080:278:C8RE3ACXX:6:2211:10900:88094|TGCT     99      chr7    55221737        60      98M     =       55222033        394
    H1080:278:C8RE3ACXX:6:2213:20942:84732|TGCT     97      chr7    55221737        60      98M     =       55222033        394

    H1080:278:C8RE3ACXX:6:2211:10900:88094|TGCT     147     chr7    55222033        60      98M     =       55221737        -394
    H1080:278:C8RE3ACXX:6:2213:20942:84732|TGCT     145     chr7    55222033        60      98M     =       55221737        -394

    In this example, location and insert size are exactly the same. Take 99 as consensus flag for first 2 reads, and
    147 for second.
    '''
    # Rank flags by number of occurrences
    count_flags = collections.Counter(i.flag for i in bam_reads).most_common()  # [(97, 1), (99, 1)]
    # List all flags with max count (will show multiple if there's a tie for the max count)
    max_flag = [i for i, j in count_flags if j == count_flags[0][1]]

    if len(max_flag) != 1:
        if max_flag == [97, 99]:
            flag = 99
        elif max_flag == [81, 83]:
            flag = 83
        elif max_flag == [145, 147]:
            flag = 147
        elif max_flag == [161, 163]:
            flag = 163
        else:
            flag = max_flag[randint(0, len(max_flag)-1)]
    else:
        flag = max_flag[0]

    return flag


def create_aligned_segment(bam_reads, sscs, sscs_qual, query_name):
    '''(list, str, list) -> pysam object
    Return pysam object with new consensus seq given list of bam reads.

    Bam file characteristics:
    1) Query name -> new 'consensus' query name (e.g. TTTG_24_58847448_24_58847416_137M10S_147M_pos_99_147)
    2) Flag -> take most common flag
    3) Reference sequence chr
    4) 1-based leftmost mapping POSition
    5) Mapping quality -> most common mapping quality
    6) Cigar string
    7) Reference sequence of mate/next read
    8) Position of mate/next read
    9) Observed template length (excluding softclips and deletions)
    10) Sequence -> consensus sequence
    11) Quality Score -> molecular consensus phred quality
    12) Additional parameters:
        - Read Group (RG) -> comb of unique read groups
        - Combination of all the reads that went into making the consensus
    '''
    # Cigar should be same across all reads
    common_cigar = bam_reads[0].cigarstring

    # Find first read with most common cigar and set as template read for SSCS
    template_index = [i.cigarstring for i in bam_reads].index(common_cigar)
    template_read = bam_reads[template_index]

    # Create bam read based on template read
    SSCS_read = pysam.AlignedSegment()
    SSCS_read.query_name = query_name
    SSCS_read.query_sequence = sscs
    SSCS_read.reference_id = template_read.reference_id
    SSCS_read.reference_start = template_read.reference_start
    SSCS_read.mapping_quality = read_mode('mapping_quality', bam_reads)  # Most common mapping quality
    SSCS_read.cigar = template_read.cigar
    SSCS_read.next_reference_id = template_read.next_reference_id
    SSCS_read.next_reference_start = template_read.next_reference_start
    SSCS_read.template_length = read_mode('template_length', bam_reads)
    SSCS_read.query_qualities = sscs_qual

    # Most common flag used unless there's a tie, then flags are ranked if its 99/97/147/145, otherwise randomly picked
    SSCS_read.flag = consensus_flag(bam_reads)

    SSCS_read.set_tag('RG', read_mode("get_tag('RG')", bam_reads))

    return SSCS_read


def reverse_seq(seq):
    '''(str) -> str
    Return reverse compliment of sequence (used for writing rev comp sequences to fastq files).

    >>> reverse_seq('TCAGCATAATT')
    'AATTATGCTGA'
    >>> reverse_seq('ACTGNN')
    'NNCAGT'
    '''
    rev_comp = ''
    nuc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for base in seq:
        rev_comp = nuc[base] + rev_comp

    return rev_comp


###############################
##           Main            ##
###############################
# if __name__ == "__main__":
#     # bamfile = pysam.AlignmentFile(
#     # '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test25/majority6/MEM-001-p03125.sscs.sort.bam',
#     # "rb")
#     #
#     # badRead_bam = pysam.AlignmentFile('/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test25/majority6/MEM-001-p03125.sscs.badReads.bam', "wb", template = bamfile)
#
#     # bamfile = pysam.AlignmentFile(
#     # '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/subsampled_bam/MEMU/MEM-001-p03125.sort.bam',
#     # "rb")
#     #
#     # badRead_bam = pysam.AlignmentFile('/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test25/majorityq7//MEM-001-p03125.badReads.bam', "wb", template = bamfile)
#
#     bamfile = pysam.AlignmentFile(
#         '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/subsampled_bam/EPIC_24norm/EPIC_edgecases.sort.bam',
#         "rb")
#
#     badRead_bam = pysam.AlignmentFile(
#         '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test25/majorityq8/EPIC_edgecases.badReads.bam',
#         "wb", template=bamfile)
#
#     # ===== Initialize dictionaries =====
#     read_dict = collections.OrderedDict()
#     tag_dict = collections.defaultdict(int)
#     pair_dict = collections.defaultdict(list)
#     csn_pair_dict = collections.defaultdict(list)
#
#     # ===== Initialize counters =====
#     unmapped = 0
#     bad_reads = 0  # secondary/supplementary reads
#     counter = 0
#     singletons = 0
#     SSCS_reads = 0
#
#     chr_data = read_bam(bamfile,
#                         pair_dict=pair_dict,
#                         read_dict=read_dict,
#                         tag_dict=tag_dict,
#                         csn_pair_dict=csn_pair_dict,
#                         badRead_bam=badRead_bam)
#
#     read_dict = chr_data[0]
#     tag_dict = chr_data[1]
#     pair_dict = chr_data[2]
#     csn_pair_dict = chr_data[3]
#
#     counter += chr_data[4]
#     unmapped += chr_data[5]
    bad_reads += chr_data[6]
