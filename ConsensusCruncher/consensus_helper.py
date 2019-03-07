#!/usr/bin/env python

###############################################################
#
#                     Consensus Helper
#
# Author: Nina Wang
# Date Created: Mar 24, 2016
###############################################################
# Function:
# Helper functions for single strand and duplex consensus making.
#
# Written for Python 3.5.1
#
# Concepts:
#   - Unique tag: identifier for grouping PCR duplicates from the same read of a strand of a molecule
#   - Consensus tag: new query name to pair consensus tags (R1 and R2 from the same strand of a molecule)
#                    Each consensus tag corresponds to 2 unique tags
#
###############################################################

##############################
#        Load Modules        #
##############################
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
    """(str) -> dict
    Return dictionary of coordinates based on bed file.
    """
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
    """(int) -> str
    Returns read number based on flag.

    Test cases:
    >>> which_read(83)
    'R1'
    >>> which_read(131)
    'R2'
    >>> which_read(177)
    'R2'
    """
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
    """(pysam.calignedsegment.AlignedSegment) -> str
    Return DNA strand of origin based on flags.
    Note: Strand is needed for the common identifier to replace read orientation and number (see sscs_qname for example)
    flag_pairings = {
                     # paired and mapped
                     99:147, 147:99, 83:163, 163:83,
                     # mapped within insert size, but wrong orientation (++, --)
                     67:131, 131:67, 115:179, 179:115,
                     # mapped uniquely, but wrong insert size
                     81:161, 161:81, 97:145, 145:97,
                     # wrong insert size and wrong orientation
                     65:129, 129:65, 113:177, 177:113
                     }
    To determine orientation of translocations (flag 65, 129, 113, 177) where read pair occurs in same direction or have
    no orientation, use coordinate to determine strand
    - pos: 1) read is first in pair (R1) AND has lower chr number than mate with higher chr number (R1 chr < R2 chr)
           2) read is second in pair (R2) AND has higher chr number than mate with lower chr number (R2 chr > R1 chr)
           3) read is first in pair (R1) AND is on same chr as mate AND start position is less than mate start
           4) read is second in pair (R2) AND is on the same chr as mate AND start position is more than mate start
           e.g. H1080:278:C8RE3ACXX:6:1209:19123:36119|CTCT	113	12	25398309	60	98M	16	7320258	98 ->
                (chr12_chr16_R1) 'pos'
                H1080:278:C8RE3ACXX:6:1209:19123:36119|CTCT	177	16	7320258	60	98M	12	25398309	98 ->
                (chr16_chr12_R2) 'pos'
    - neg: Opposite cases of pos (R1 chr > R2 chr | R2 chr < R1 chr)
           e.g. H1080:278:C8RE3ACXX:6:1307:20616:55254|CTCT	177	12	25398309	60	98M	16	7320258	98 ->
                (chr16_chr12_R1) 'neg'
                H1080:278:C8RE3ACXX:6:1307:20616:55254|CTCT	113	16	7320258	60	98M	12	25398309	98 ->
                (chr12_chr16_R2) 'neg'
    Exceptions:
    81/161 - reads from two strands of a molecule encoded as 81/161/81/161 when typically different flags given e.g.
             99/147/83/163 between two strands. In addition, duplex pairing of these reads also differ (ignore edge
             cases for now as there are very few occurrences)
    normal:
        for cases where flags of read/mate and their reverse complement are encoded in
          the same direction or share the same flags
          (e.g. flags 65/129 -> +/+/+/+ OR 81/161/81/161 [scenarios like this cannot be distinguished by flag info such
          as orientation and strand number alone. Genome coordinate can be included to help differentiate reads])
    Example Test cases:
    Flag = 147 -> 'pos'
    Flag = 67 -> 'pos'
    Flag = 131 -> 'pos'
    Flag = 81 -> 'neg'
    """
    # Flags indicating strand direction
    pos = [99, 147, 67, 131]
    neg = [83, 163, 115, 179]
    no_ori = [65, 129, 113, 177, 81, 161, 97, 145]  # direction not defined

    if read.flag in pos:
        strand = 'pos'
    elif read.flag in neg:
        strand = 'neg'
    elif read.flag in no_ori:
        # Determine orientation of flags with no defined direction using order of chr coor
        if (read.reference_id < read.next_reference_id and which_read(read.flag) == 'R1') or \
                (read.reference_id > read.next_reference_id and which_read(read.flag) == 'R2') or \
                (read.reference_id == read.next_reference_id and which_read(read.flag) == 'R1' and
                         read.reference_start < read.next_reference_start) or \
                (read.reference_id == read.next_reference_id and which_read(read.flag) == 'R2' and
                         read.reference_start > read.next_reference_start):
            strand = 'pos'
        else:
            strand = 'neg'
    else:
        # Only uniquely mapped reads (with flags indicated above) should be retained, as 'bad reads' were filtered out
        # in a previous step
        print('STRAND ERROR')
        print(read.flag)
        strand = None

    return strand


def cigar_order(read, mate):
    """(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> str
    Return ordered cigar string from paired reads based on strand and read number.
    * Note: order does not correspond to read and mate as cigars were extracted from read pair prior to assignment of
    individual read identifiers (as mate cigar is not a parameter of each read, we would have to track both cigars
    individually and order them for each read separately if we wanted cigars to be in the order of mate_pair. It just
    makes more intuitive sense to track them by strand and read number since that's what we'll be using for the new
    query name and that's the part end-users see anyways)

    0 = Match/mismatch
    1 = Insertion
    2 = Deletion
    4 = Soft clip
    5 = Hard clip

    Pos strand, R1 cigar first
    Neg strand, R2 cigar first
    => want cigar in same order for Duplex consensus sequence tag identifier

    Examples:
        (+) [99]  137M10S    => 137M10S_147M
            [147] 147M

        (-) [83]  147M   => 137M10S_147M
            [163] 137M10S
    """
    ori_strand = which_strand(read)
    read_num = which_read(read.flag)

    if (ori_strand == 'pos' and read_num == 'R1') or (ori_strand == 'neg' and read_num == 'R2'):
        cigar = '{}_{}'.format(read.cigarstring,
                               mate.cigarstring)
    else:
        cigar = '{}_{}'.format(mate.cigarstring,
                               read.cigarstring)

    return cigar


def sscs_qname(read, mate, barcode, cigar):
    """(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> str
    Return new query name for consensus sequences:
    [Barcode]_[Read Chr]_[Read Start]_[Mate Chr]_[Mate Start]_[Read Cigar String]_[Mate Cigar String]_[Strand]_[Absolute insert size]:[Family Size]

    * Since multiple reads go into making a consensus, a new query name is needed as an identifier for consensus read
    pairs * (Read pairs share the same query name to indicate they're mates)

    Note: coordinates are ordered from low -> high and cigars ordered by read number and strand, so read pairs will
          share common identifiers
    WARNING: Pysam converts coordinates to be 0-based, query name coordinates don't match coordinate seen in read

    ---
    original query name -> H1080:278:C8RE3ACXX:6:1308:18882:18072|TTTG
    new query name -> TTTG_24_58847448_24_58847416_137M10S_147M_pos_148

    Examples:
    (+)                                                 [Flag]
    TTTG_24_58847416_24_58847448_137M10S_147M_pos_fwd_R1 [99] --> TTTG_24_58847416_24_58847448_137M10S_147M_pos_148
    TTTG_24_58847448_24_58847416_137M10S_147M_pos_rev_R2 [147]

    (-)
    TGTT_24_58847448_24_58847416_137M10S_147M_neg_rev_R1 [83] --> TGTT_24_58847416_24_58847448_137M10S_147M_neg_148
    TGTT_24_58847416_24_58847448_137M10S_147M_neg_fwd_R2 [163]

    Special case (mate and complementary reads are all in the same direction)
    - Use coordinate and flags to differentiate between strand (see which_strand fx for details)
    """
    read_chr = read.reference_id
    mate_chr = mate.reference_id
    read_coor = read.reference_start
    mate_coor = mate.reference_start

    if (read_chr == mate_chr and int(read_coor) > int(mate_coor)) or \
            (int(read_chr) > int(mate_chr)):
        read_chr = mate.reference_id
        mate_chr = read.reference_id
        read_coor = mate.reference_start
        mate_coor = read.reference_start

    strand = which_strand(read)
    query_tag = '{}_{}_{}_{}_{}_{}_{}_{}'.format(barcode,
                                                 read_chr,
                                                 read_coor,
                                                 mate_chr,
                                                 mate_coor,
                                                 cigar,
                                                 strand,
                                                 abs(read.template_length))

    return query_tag


def unique_tag(read, barcode, cigar):
    """(pysam.calignedsegment.AlignedSegment, str, str) -> str
    Return unique identifier tag for one read of a strand of a molecule.

    Tag uses following characteristics to group reads belonging to the same strand of an individual molecule (PCR dupes):
    [Barcode]_[Read Chr]_[Read Start]_[Mate Chr]_[Mate Start]_[Cigar String]_[Orientation]_[ReadNum]
    e.g. TTTG_24_58847416_24_58847448_137M10S_147M_fwd_R1

    Notes:
        - barcode always in order of R1/R2 for each read pair as determined from tag_to_header.py extraction of tags
          from fastq files
        - paired reads have ordered cigar strings of read and mate (see cigar_order fx for details) - important to have
        both for easy duplex tag search (if we didn't track both, you'd have to look for the corresponding mate cigar
        each time)
        - barcode location differs depending on whether you're making SSCS (end of header in uncollapsed) vs
          DCS (start of header in SSCS)
        - chr of R1 and R2 included as they might differ in case of translocation
        - orientation and read number included to differentiate palindromic barcodes

    Example:
       R1 --->   <--- R2
    (+) TT-----------GT

    (-) TT-----------GT
       R2 --->   <--- R1

    R1 of (+) -> TTTG_24_58847416_24_58847448_137M10S_147M_fwd_R1
    NB500964:12:HTTG2BGXX:4:22601:26270:1144|TTTG	99	24	58847416	17	137M10S	24	58847448	137

    R2 of (+) -> TTTG_24_58847448_24_58847416_137M10S_147M_rev_R2
    NB500964:12:HTTG2BGXX:4:22601:26270:1144|TTTG	147	24	58847448	17	147M	24	58847416	147

    R1 of (-) -> TGTT_24_58847448_24_58847416_137M10S_147M_rev_R1

    R2 of (-) -> TGTT_24_58847416_24_58847448_137M10S_147M_fwd_R2
    """
    orientation = 'fwd'
    if read.is_reverse:
        orientation = 'rev'

    readNum = which_read(read.flag)

    # Unique identifier for strand of individual molecules
    tag = '{}_{}_{}_{}_{}_{}_{}_{}'.format(barcode,  # mol barcode
                                           read.reference_id,  # chr
                                           read.reference_start,  # start (0-based)
                                           read.next_reference_id,  # mate chr
                                           read.next_reference_start,  # mate start
                                           cigar,
                                           orientation,  # strand direction
                                           readNum
                                           )
    return tag


def read_bam(bamfile, pair_dict, read_dict, csn_pair_dict, tag_dict, badRead_bam, duplex,
             read_chr=None, read_start=None, read_end=None, barcode_delim=None):
    """(bamfile, dict, dict, dict, dict, bamfile, bool, str, int, int) ->
    dict, dict, dict, dict, int, int, int

    === Input ===
    - bamfile (pysam.AlignmentFile object): uncollapsed BAM file

    - pair_dict: dictionary of paired reads based on query name to process data in pairs

    - read_dict: dictionary of bamfile reads grouped by unique molecular tags

    - csn_pair_dict: dictionary of paired tags sharing the same consensus tag to track pairing

    - tag_dict: integer dictionary indicating number of reads in each read family
                 {read_tag: 2, ..etc}

    - badRead_bam (pysam.AlignmentFile object): BAM file of "bad" reads (unmapped, multiple mapping)

    -- Optional --
    # For large bamfiles that are split into regions
    - read_chr (str): chromosome region to fetch reads
    - read_start (int): starting position to fetch reads
    - read_end (int): stopping position to fetch reads

    # For duplex consensus making
    - duplex: any string or bool [that is not None] specifying duplex consensus making [e.g. TRUE], necessary for
              parsing barcode as query name for Uncollapsed and SSCS differ

    # For bams with barcodes extracted by other software and placed into read name with different delimiters
    - barcode_delim (str): sequence before barcode (e.g. '|' for 'HWI-D00331:196:C900FANXX:7:1110:14056:43945|TTTT')

    === Output ===
    1) read_dict: dictionary of bamfile reads grouped by unique molecular tags
                  Example: {read_tag: [<pysam.calignedsegment.AlignedSegment>, <pysam.calignedsegment.AlignedSegment>]}
                  - Key: [Barcode]_[Read Chr]_[Read Start]_[Mate Chr]_[Mate Start]_[Cigar String]_[Strand]_[Orientation]_[ReadNum]
                  - Value: List of reads (pysam object)

    2) tag_dict: integer dictionary indicating number of reads in each read family
                 {read_tag: 2, ..etc}

    3) pair_dict: dictionary of paired reads based on query name to process data in pairs
                 (note: this is a tmp dict as values are removed from dict once pair assigned to other dicts, this is
                 important for retaining data from translocations or reads crossing bam division regions)
                 Example: {query name: [read 1, read 2]}

    4) csn_pair_dict: dictionary of paired tags sharing the same consensus tag to track pairing (paired reads share the
                     same query name/header)
                     Example: {consensus_tag: [R1_tag, R2_tag]}

    5) counter: total number of reads

    6) unmapped: unmapped reads

    7) multiple_mapping: number of reads that not properly mapped
                         - secondary reads: same sequence aligns to multiple locations
                         - supplementary reads: multiple parts of sequence align to multiple locations
    """
    # Fetch data given genome coordinates
    if read_chr is None:
        bamLines = bamfile.fetch(until_eof=True)
    else:
        bamLines = bamfile.fetch(read_chr, read_start, read_end)

    # Initialize counters
    unmapped = 0
    unmapped_mate = 0
    multiple_mapping = 0  # secondary/supplementary reads
    counter = 0
    bad_spacer = 0

    for line in bamLines:
        # Parse out reads that don't fall within region
        if read_chr is not None:
            # pysam fetch will retrieve reads that fall outside region due to pairing (we filter out to prevent double
            # counting as we'll be fetching those reads again when we iterate through the next region)
            if line.reference_start < read_start or line.reference_start > read_end:
                continue

        counter += 1

        ######################
        #    Filter Reads    #
        ######################
        # === 1) FILTER OUT UNMAPPED / MULTIPLE MAPPING READS ===
        mate_unmapped = [73, 89, 121, 153, 185, 137]
        badRead = True

        # Check if delimiter is found in read
        if barcode_delim is not None and barcode_delim not in line.qname:
            bad_spacer += 1
        elif line.is_unmapped:
            unmapped += 1
            counter -= 1
        elif line.flag in mate_unmapped:
            unmapped_mate += 1
        elif line.is_secondary:
            multiple_mapping += 1
        elif line.is_supplementary:
            multiple_mapping += 1
        else:
            badRead = False

        # Write bad reads to file
        if badRead and badRead_bam is not None:
            badRead_bam.write(line)
        else:
            pair_dict[line.qname].append(line)

            ######################
            #      Unique ID     #
            ######################
            # === 2) ASSIGN UNIQUE IDENTIFIER TO READ PAIRS ===
            if len(pair_dict[line.qname]) == 2:
                read = pair_dict[line.qname][0]
                mate = pair_dict[line.qname][1]
                # === Create consensus identifier ===
                # Extract molecular barcode, barcodes in diff position for SSCS vs DCS generation
                if duplex == None or duplex == False:
                    if barcode_delim is None:
                        # SSCS query name: H1080:278:C8RE3ACXX:6:1308:18882:18072|CACT
                        barcode = read.qname.split("|")[1]
                    else:
                        barcode = read.qname.split(barcode_delim)[1]
                else:
                    # DCS query name: CCTG_12_25398000_12_25398118_neg:5
                    barcode = read.qname.split("_")[0]

                # Consensus_tag cigar (ordered by strand and read)
                cigar = cigar_order(read, mate)
                # Assign consensus tag as new query name for paired consensus reads
                consensus_tag = sscs_qname(read, mate, barcode, cigar)

                for i in range(2):
                    read_i = pair_dict[line.qname][i]
                    # Molecular identifier for grouping reads belonging to the same read of a strand of a molecule
                    tag = unique_tag(read_i, barcode, cigar)

                    ######################
                    #   Assign to Dict   #
                    ######################
                    # === 3) ADD READ PAIRS TO DICTIONARIES ===
                    if tag not in read_dict and tag not in tag_dict:
                        read_dict[tag] = [read_i]
                        tag_dict[tag] += 1

                        # Group paired unique tags using consensus tag
                        if consensus_tag not in csn_pair_dict:
                            csn_pair_dict[consensus_tag] = [tag]
                        elif len(csn_pair_dict[consensus_tag]) == 2:
                            # Honestly this shouldn't happen anymore with these identifiers
                            print("Consensus tag NOT UNIQUE -> multiple tags (4) share same consensus tag [due to poor "
                                  "strand differentiation as a result of identifiers lacking complexity]")
                            print(consensus_tag)
                            print(tag)
                            print(read_i)
                            print(csn_pair_dict[consensus_tag])
                            print(read_dict[csn_pair_dict[consensus_tag][0]][0])
                            print(read_dict[csn_pair_dict[consensus_tag][1]][0])

                            # Manual inspection should be done on these reads
                        else:
                            csn_pair_dict[consensus_tag].append(tag)
                    elif tag in tag_dict and read not in read_dict[tag]:
                        # Append reads sharing the same unique tag together (PCR dupes)
                        read_dict[tag].append(read_i)
                        tag_dict[tag] += 1
                    else:
                        # Data fetch error - line read twice (if its found in tag_dict and read_dict)
                        print('Pair already written: line read twice - check to see if read overlapping / near cytoband'
                              ' region (point of data division)')

                # remove read pair qname from pair_dict once reads added to read_dict
                pair_dict.pop(line.qname)

    return read_dict, tag_dict, pair_dict, csn_pair_dict, counter, unmapped_mate, multiple_mapping, bad_spacer


def read_mode(field, bam_reads):
    """(str, lst) -> str
    Return mode (most common occurrence) of a specified field

    Field e.g. cigarstring, flag, mapping quality, template_length
    """
    field = 'i.{}'.format(field)
    # Rank by number of occurrences
    field_lst = collections.Counter(eval(field) for i in bam_reads).most_common()
    # Take max occurrences
    common_field_lst = [i for i, j in field_lst if j == field_lst[0][1]]
    # Randomly select max if there's multiple
    common_field = common_field_lst[randint(0, len(common_field_lst)-1)]

    return common_field


def consensus_flag(bam_reads):
    """(list) -> str
    Return consensus flag given list of reads from the same family.

    If multiple flags are present within reads from the same molecule and a max can't be determined, prioritize flags
    that indicate proper mapping/pairing (99, 147, 83, 163) and within insert size.
    e.g.
    H1080:278:C8RE3ACXX:6:2211:10900:88094|TGCT     99      chr7    55221737        60      98M     =       55222033        394
    H1080:278:C8RE3ACXX:6:2213:20942:84732|TGCT     97      chr7    55221737        60      98M     =       55222033        394

    H1080:278:C8RE3ACXX:6:2211:10900:88094|TGCT     147     chr7    55222033        60      98M     =       55221737        -394
    H1080:278:C8RE3ACXX:6:2213:20942:84732|TGCT     145     chr7    55222033        60      98M     =       55221737        -394

    In this example, location and insert size are exactly the same. Take 99 as consensus flag for first 2 reads, and
    147 for second.
    """
    # Rank flags by number of occurrences
    count_flags = collections.Counter(i.flag for i in bam_reads).most_common()  # [(97, 1), (99, 1)]
    # List all flags with max count (will show multiple if there's a tie for the max count)
    max_flag = [i for i, j in count_flags if j == count_flags[0][1]]

    if len(max_flag) != 1:
        if 99 in max_flag:
            flag = 99
        elif 83 in max_flag:
            flag = 83
        elif 147 in max_flag:
            flag = 147
        elif 163 in max_flag:
            flag = 163
        else:
            flag = max_flag[randint(0, len(max_flag)-1)]  # If flag not properly paired/mapped, randomly select from max
    else:
        flag = max_flag[0]

    return flag


def create_aligned_segment(bam_reads, sscs, sscs_qual, query_name):
    """(list, str, list, list, str) -> pysam object
    Return consensus read representing list of reads from the same molecule.

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
        - Read Group (RG)
        - Proportion score (ps) -> we calculated in consensus_maker
            ->  Tags starting with ‘X’, ‘Y’ or ‘Z’ and tags containing lowercase letters in either position are reserved
                for local use and will not be formally defined in any future version of these specifications.
    """
    # Use first read in list as template (all reads should share same cigar, template length, and coor)
    template_read = bam_reads[0]

    # Create consensus read based on template read
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

    # Most common flag used unless there's a tie, then flags are ranked if its 99/83/147/163, otherwise randomly picked
    SSCS_read.flag = consensus_flag(bam_reads)

    # Optional fields
    try:
        SSCS_read.set_tag('RG', read_mode("get_tag('RG')", bam_reads))
    except:
        pass

    return SSCS_read


def reverse_seq(seq):
    """(str) -> str
    Return reverse complement of sequence (used for writing rev comp sequences to fastq files).

    >>> reverse_seq('TCAGCATAATT')
    'AATTATGCTGA'
    >>> reverse_seq('ACTGNN')
    'NNCAGT'
    """
    rev_comp = ''
    nuc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for base in seq:
        rev_comp = nuc[base] + rev_comp

    return rev_comp


def duplex_tag(tag):
    """(str) -> str
    Return tag for duplex read.

    Things to be changed in tag to find its complementary tag:
    1) barcode: molecular identifiers get swapped (e.g. 2 based identifiers on each side of DNA fragment)
               (+) 5' AT-------CG  3' -> ATGC
               (-)    AT-------CG     <- GCAT
    2) read: R1 -> R2

    Note: don't need to swap cigar strings as they are already ordered by strand (pos R1 correspond to neg R2)
    ** Barcode lists may contain barcodes of different lengths, so R1 and R2 barcodes are 
    separated by '.'**

    Test cases:
    >>> duplex_tag('GTCT_1_1507809_7_55224319_98M_98M_fwd_R1')
    'CTGT_1_1507809_7_55224319_98M_98M_neg_fwd_R2'
    >>> duplex_tag('GTCT_7_55224319_1_1507809_98M_98M_rev_R2')
    'CTGT_7_55224319_1_1507809_98M_98M_neg_rev_R1'
    >>> duplex_tag('CTGT_1_1507809_7_55224319_98M_98M_fwd_R2')
    'GTCT_1_1507809_7_55224319_98M_98M_pos_fwd_R1'
    >>> duplex_tag('CTGT_7_55224319_1_1507809_98M_98M_rev_R1')
    'GTCT_7_55224319_1_1507809_98M_98M_pos_rev_R2'
    """
    split_tag = tag.split('_')
    # 1) Barcode needs to be swapped
    barcode = split_tag[0]
    # Separate R1 and R2 barcodes with '.' separator
    if re.search('\.', barcode) is not None:
        split_index = barcode.index('.')
        split_tag[0] = barcode[split_index+1:] + '.' + barcode[:split_index]
    else:
        barcode_bases = int(len(barcode) / 2)  # number of barcode bases, avoids complications if num bases change
        # duplex barcode is the reverse (e.g. AT|GC -> GC|AT [dup])
        split_tag[0] = barcode[barcode_bases:] + barcode[:barcode_bases]
        		
    # 2) Opposite read number in duplex
    read_num = split_tag[8]
    if read_num == 'R1':
        split_tag[8] = 'R2'
    else:
        split_tag[8] = 'R1'

    return '_'.join(split_tag)
