#!/usr/bin/env python3

###############################################################
#
#           Duplex Consensus Sequence (DCS) Generator
#
# Author: Nina Wang
# Date Created: Mar 24, 2016
###############################################################
# Function:
# To generate duplex/double-strand consensus sequences for molecule based error suppression.
# - Consensus sequence from bases with Phred quality >= 30
# - Consensus quality score from addition of quality scores (i.e. product of error probabilities)
#
# Written for Python 3.5.1
#
# Usage:
# Python3 DCS_maker.py [--infile INFILE] [--outfile OUTFILE] [--bedfile BEDFILE]
#
# Arguments:
# --infile INFILE     input BAM file
# --outfile OUTFILE   output BAM file
# --bedfile BEDFILE   Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt -
#                     See bed_separator.R for making your own bed file based on a target panel / specific coordinates)
#
# Inputs:
# 1. A position-sorted BAM file containing paired-end reads with SSCS consensus identifier in the header/query name
# 2. A BED file containing coordinates subdividing the entire ref genome for more manageable data processing
#
# Outputs:
# 1. A BAM file containing paired double stranded consensus sequences - "dcs.bam"
# 2. A SSCS singleton BAM file containing SSCSs without reads from the complementary strand - "sscs.singleton.bam"
# 3. A text file containing summary statistics (Total SSCS reads, Unmmaped SSCS reads, Secondary/Supplementary SSCS
#    reads, DCS reads, and SSCS singletons) - "stats.txt" (Stats pended to same stats file as SSCS)
#
# Concepts:
#    - Read family: reads that share the same molecular barcode, chr, and start
#                   coordinates for Read1 and Read2
#    - SSCS Singleton: a SSCS read without its complementary read
#
###############################################################

##############################
#        Load Modules        #
##############################
import pysam  # Need to install
import collections
import re
import array
from random import randint
from argparse import ArgumentParser
import math

from consensus_helper import *


###############################
#       Helper Functions      #
###############################
def dcs_consensus_tag(tag, ds):
    """(str, str) -> str
    Return consensus tag for duplex reads.

    Strand removed and family size from both SSCS included in order of pos_neg strand.

    Test cases:
    >>> dcs_consensus_tag('TTCA_7_55259315_7_55259454_98M_98M_neg:3', 'CATT_7_55259315_7_55259454_98M_98M_pos:6')
    'CATT_TTCA_7_55259315_7_55259454_98M_98M:6_3'

    >>> dcs_consensus_tag('CTTC_23_74804535_23_74804611_98M_98M_pos:2', 'TCCT_23_74804535_23_74804611_98M_98M_neg:3')
    'CTTC_TCCT_23_74804535_23_74804611_98M_98M:2_3'

    >>> dcs_consensus_tag('TTTC_7_140477735_7_140477790_98M_98M_neg:3', 'TCTT_7_140477735_7_140477790_98M_98M_pos:2')
    'TCTT_TTTC_7_140477735_7_140477790_98M_98M:2_3'
    """
    barcode = tag.split('_')[0]
    duplex_barcode = ds.split('_')[0]
    tag_coor = tag.split('_', 1)[1].rsplit('_', 1)[0]
    tag_fam_size = tag.split(':')[1]
    ds_fam_size = ds.split(':')[1]

    # Order tag barcodes and family size based on strand (pos then negative)
    if 'pos' in tag:
        dcs_query_name = "{}_{}_{}:{}_{}".format(barcode,
                                                 duplex_barcode,
                                                 tag_coor,
                                                 tag_fam_size,
                                                 ds_fam_size)
    else:
        dcs_query_name = "{}_{}_{}:{}_{}".format(duplex_barcode,
                                                 barcode,
                                                 tag_coor,
                                                 ds_fam_size,
                                                 tag_fam_size)

    return dcs_query_name


def duplex_consensus(read1, read2):
    """(pysam.calignedsegment.AlignedSegment, pysam.calignedsegment.AlignedSegment) -> pysam.calignedsegment.AlignedSegment

    Return consensus of complementary reads with N for inconsistent bases.
    """
    consensus_seq = ''
    consensus_qual = []

    for i in range(read1.query_length):
        # Check to see if base at position i is the same
        if read1.query_sequence[i] == read2.query_sequence[i]:
            consensus_seq += read1.query_sequence[i]
            mol_qual = sum([read1.query_qualities[i], read2.query_qualities[i]])
            # Set to max quality score if sum of qualities is greater than the threshold (Q60) imposed by genomic tools
            if mol_qual > 60:
                consensus_qual += [60]
            else:
                consensus_qual += [mol_qual]
        else:
            consensus_seq += 'N'
            consensus_qual += [0]

    return consensus_seq, consensus_qual


###############################
#        Main Function        #
###############################

def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--infile", action = "store", dest="infile", help="Input BAM file", required=True)
    parser.add_argument("--outfile", action = "store", dest="outfile", help="Output BAM file", required=True)
    parser.add_argument("--bedfile", action="store", dest="bedfile",
                        help="Bedfile containing coordinates to subdivide the BAM file (Recommendation: cytoband.txt - \
                        See bed_separator.R for making your own bed file based on a target panel/specific coordinates)",
                        required=False)
    args = parser.parse_args()

    ######################
    #       SETUP        #
    ######################
    start_time = time.time()
    # ===== Initialize input and output bam files =====
    args.infile = str(args.infile)
    args.outfile = str(args.outfile)

    sscs_bam = pysam.AlignmentFile(args.infile, "rb")
    dcs_bam = pysam.AlignmentFile(args.outfile, "wb", template=sscs_bam)
    
    if re.search('dcs\.sc', args.outfile) is not None:
        sscs_singleton_bam = pysam.AlignmentFile('{}.sscs.sc.singleton.bam'.format(args.outfile.split('.dcs.sc')[0]),
                                             "wb", template=sscs_bam)
        dcs_header = "DCS - Singleton Correction"
        sc_header = " SC"
    else:
        sscs_singleton_bam = pysam.AlignmentFile('{}.sscs.singleton.bam'.format(args.outfile.split('.dcs')[0]),
                                             "wb", template=sscs_bam)
        dcs_header = "DCS"
        sc_header = ""

    stats = open('{}.stats.txt'.format(args.outfile.split('.dcs')[0]), 'a')
    time_tracker = open('{}.time_tracker.txt'.format(args.outfile.split('.dcs')[0]), 'a')

    # ===== Initialize dictionaries and counters=====
    read_dict = collections.OrderedDict()
    tag_dict = collections.defaultdict(int)
    pair_dict = collections.defaultdict(list)
    csn_pair_dict = collections.defaultdict(list)

    unmapped = 0
    unmapped_mate = 0
    multiple_mapping = 0  # Secondary/supplementary reads
    counter = 0
    sscs_singletons = 0  # Single strand consensus sequences without a complementary strand
    multiple_mappings = 0

    duplex_count = 0
    duplex_dict = collections.defaultdict(int)

    #######################
    #   SPLIT BY REGION   #
    #######################
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
            read_chr = x.split('_', 1)[0]
            read_start = division_coor[x][0]
            read_end = division_coor[x][1]

        chr_data = read_bam(sscs_bam,
                            pair_dict=pair_dict,
                            read_dict=read_dict,
                            csn_pair_dict=csn_pair_dict,
                            tag_dict=tag_dict,
                            badRead_bam=None,
                            duplex=True,
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
        multiple_mapping += chr_data[6]

        ######################
        #     CONSENSUS      #
        ######################
        # ===== Create consenus seq for reads =====
        for readPair in list(csn_pair_dict.keys()):
            for tag in csn_pair_dict[readPair]:
                # Determine tag of duplex read
                ds = duplex_tag(tag)

                # === Group duplex read pairs and create consensus ===
                # Check presence of duplex pair
                if ds not in duplex_dict.keys():
                    if tag in tag_dict and ds in tag_dict:
                        duplex_count += 1

                        # consensus seq
                        consensus_seq, consensus_qual = duplex_consensus(read_dict[tag][0], read_dict[ds][0])

                        # consensus duplex tag
                        dcs_query_name = dcs_consensus_tag(read_dict[tag][0].qname, read_dict[ds][0].qname)  # New query name containing both barcodes

                        dcs_read = create_aligned_segment([read_dict[tag][0], read_dict[ds][0]], consensus_seq,
                                                          consensus_qual, dcs_query_name)

                        # add duplex tag to dictionary to prevent making a duplex for the same sequences twice
                        duplex_dict[tag] += 1

                        dcs_bam.write(dcs_read)

                    else:
                        sscs_singleton_bam.write(read_dict[tag][0])
                        sscs_singletons += 1

                    # Remove read from dictionary after writing
                    del read_dict[tag]

            # Remove key from dictionary after writing
            del csn_pair_dict[readPair]

    ######################
    #       SUMMARY      #
    ######################
    summary_stats = '''# === {} ===
SSCS{} - Total reads: {}
SSCS{} - Unmapped reads: {}
SSCS{} - Secondary/Supplementary reads: {}
DCS{} reads: {}
SSCS{} singletons: {} \n'''.format(dcs_header, sc_header, counter, sc_header, unmapped, sc_header, multiple_mappings,
                                   sc_header, duplex_count, sc_header, sscs_singletons)
    stats.write(summary_stats)
    print(summary_stats)

    # Output total DCS time
    time_tracker.write('DCS: ')
    time_tracker.write(str((time.time() - start_time)/60) + '\n')

    # Close files
    time_tracker.close()
    stats.close()
    dcs_bam.close()
    sscs_singleton_bam.close()

    return duplex_dict


###############################
#            Main             #
###############################
if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
