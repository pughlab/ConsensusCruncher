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


def which_strand(flag):
    '''(int) -> str
    Return DNA strand based on flags.

    Test cases:
    >>> which_strand(147)
    'pos'
    >>> which_strand(67)
    'pos'
    >>> which_strand(131)
    'pos'
    >>> which_strand(81)
    'neg'
    '''
    pos = [99, 147, 67, 131, 97, 145, 65, 129]  # note 65/129 direction not defined
    neg = [83, 163, 115, 179, 81, 161, 113, 177]

    if flag in pos:
        strand = 'pos'
    elif flag in neg:
        strand = 'neg'
    else:
        # Only uniquely mapped reads (with flags indicated above) should be retained as 'bad reads' are filtered out
        print('STRAND ERROR')
        print(flag)
        strand = None

    return strand


def cigar_order(read, read_pair):
    '''(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> str
    Return ordered cigar string from R1 and R2

    Pos strand, R1 cigar first
    Neg strand, R2 cigar first
    => want cigar in same order for Duplex consensus formation

    Examples:
        (+) [99]  137M10S    => 137M10S_147M
            [147] 147M

        (-) [83]  147M   => 137M10S_147M
            [163] 137M10S
    '''
    ori_strand = which_strand(read.flag)
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

    strand = which_strand(read.flag)

    orientation = 'fwd'
    if read.is_reverse:
        orientation = 'rev'

    readNum = which_read(read.flag)

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
    Return new queryname for consensus sequence: barcode_chr_start_chr_end_strand_flags

    * Since multiple reads go into making a consensus, a new unique identifier
    is required to match up the read with its pair *

    Reason why we want both the flag and strand:
    - flag: to differentiate between reads with complimentary read in same orientation (e.g. 65 / 129)
    - strand: needed to make duplex consensus sequences as reads are grouped based on strand and not specific flags

    Note: coordinates are ordered from low -> high

    Examples:
    (+)                                                 [Flag]
    TTTG_24_58847416_24_58847448_137M10S_147M_pos_fwd_R1 [99] --> TTTG_24_58847448_24_58847416_137M10S_147M_pos_99_147
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
    'TGTT_24_58847448_24_58847416_137M10S_147M_neg_83_163'
    >>> sscs_qname('TGTT_24_58847416_24_58847448_137M10S_147M_neg_fwd_R2', 163)
    'TGTT_24_58847448_24_58847416_137M10S_147M_neg_83_163'
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
    # Need pos/neg to differentiate between reads with same flag from diff strands (e.g. 113, 177)
    # with smaller flag ordered first, to help differentiate between reads
    # (e.g. read and its duplex might have same coordinate and strand direction,
    # after ordering coordinates from smallest -> biggest)
    if flag < flag_pairings[flag]:
        new_tag = '{}_{}_{}'.format(new_tag, flag, flag_pairings[flag])
    else:
        new_tag = '{}_{}_{}'.format(new_tag, flag_pairings[flag], flag)

    return new_tag


def read_bam(bamfile, pair_dict, read_dict, tag_dict, badRead_bam, read_chr = None,
             read_start = None, read_end = None, duplex = None):
    '''(bamfile object, dict, dict, dict, str, int, int, str) ->
    dict, dict, dict, dict, dict, int, int, int, int

    === Input ===
    - bamfile: pysam.AlignmentFile object of bamfile

    - pair_dict: dictionary of paired tags (retains data from translocations or reads crossing cytobands to preserve
                 pairing as bamfiles are divided into sections)
                 {query name: [read 1, read 2]}

    - read_dict: dictionary of reads sharing a common tag (aka from the same unique molecule)
                 {read_tag: [<pysam.calignedsegment.AlignedSegment>,
                  <pysam.calignedsegment.AlignedSegment>, ..etc.]}

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
    if read_chr == None:
        bamLines = bamfile.fetch(until_eof = True)
    else:
        bamLines = bamfile.fetch(read_chr, read_start, read_end)

    unmapped = 0
    unmapped_flag = 0
    bad_reads = 0 # secondary/supplementary reads
    poor_mapq = 0
    counter = 0

    for line in bamLines:
        counter += 1

        # === Filter out 'bad' reads ====
        # Unmapped flags
        bad_flags = [73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141]
        badRead = True
        if line.is_unmapped:  # flag != 0
            unmapped += 1
        elif line.is_secondary:
            bad_reads += 1
        elif line.is_supplementary:
            bad_reads += 1
        elif line.flag in bad_flags:
            unmapped_flag += 1
        elif line.mapping_quality <= 5:
            poor_mapq += 1
        else:
            badRead = False

        # Write bad reads to file
        if badRead:
            badRead_bam.write(line)
        else:
            # === Add reads to dictionary as pairs ===
            if line.qname in line.qname and line in pair_dict[line.qname]:
                print('Read already read once!')
                print(read_chr, read_start, read_end)
                print(line.qname)
                print(line)
                print(pair_dict[line.qname])
                counter -= 1
            else:
                pair_dict[line.qname].append(line)

                # Create tag once reads are paired
                if len(pair_dict[line.qname]) == 2:
                    cigar = cigar_order(pair_dict[line.qname][0], pair_dict[line.qname][1])

                    for read in pair_dict[line.qname]:
                        # Barcodes extracted from diff position for duplex consensus formation
                        if duplex == None or duplex == False:
                            # SSCS query name: H1080:278:C8RE3ACXX:6:1308:18882:18072|CACT
                            barcode = read.qname.split("|")[1]
                        else:
                            # DCS query name: CCTG_12_25398000_12_25398118_neg_83_163:5
                            barcode = read.qname.split("_")[0]

                        # Unique identifier for strand of individual molecules
                        tag = unique_tag(read, cigar, barcode)

                        # Create consensus tag to be used as new query name (same query name between 2 reads indicate
                        # they're mate pairs)
                        consensus_tag = sscs_qname(tag, int(read.flag))

                        if tag not in read_dict and tag not in tag_dict:
                            read_dict[tag] = [read]
                        elif tag in tag_dict and read not in read_dict[tag]:
                            read_dict[tag].append(read)
                        else:
                            # If tag found in tag_dict, but not read_dict means line was previously read and written
                            # to file (hence removal from read_dict but presence in tag_dict)
                            print('Pair already written: line read twice - double check to see if its overlapping / near cytoband region (point of data division)')
                            print(read_chr, read_start, read_end)
                            print(tag)
                            print(read)
                            print(tag_dict[tag])
                            print(read_dict[tag][0])
                            print(consensus_tag)
                            print(tag in tag_dict)
                            print(pair_dict[consensus_tag])
                            counter -= 1
                            continue

                        tag_dict[tag] += 1

                elif len(pair_dict[line.qname]) == 3:
                    print('ERROR: triplet reads with same qname????')
                    print(line)

                    # remove read pair qname from pair_dict once reads added to read_dict as pairs
                    pair_dict.pop(line.qname)

    return read_dict, tag_dict, pair_dict, counter, unmapped, unmapped_flag, \
           bad_reads, poor_mapq


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
