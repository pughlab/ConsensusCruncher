#!/usr/bin/env python3

###############################################################
#
#                       Extract Barcodes
#
###############################################################
# Function:
# To isolate duplex barcodes from paired-end sequence reads and store in FASTQ headers after removal of spacer regions.
#
# Written for Python 3.5.1
#
# USAGE:
# python3 extract_barcodes.py [--read1 READ1] [--read2 READ2] [--outfile OUTFILE] [--blen BARCODELEN] [--slen SPACERLEN]
#                            [--sfilt SPACERFILT]
#
# Arguments:
# --read1 READ1         Input FASTQ file for Read 1 (unzipped)
# --read2 READ2         Input FASTQ file for Read 2 (unzipped)
# --outfile OUTFILE     Output FASTQ files for Read 1 and Read 2 using given filename
# --blen BARCODELEN     Barcode length
# --slen SPACERLEN      Spacer length (This region is removed and only the barcode will be added to the header)
# --sfilt SPACERFILT    Filter that excludes reads without the specified base(s) in the spacer region
# --ilen INDENTLEN      Indent length (This region is removed and only the barcode will be added to the header)
# --ifilt INDENTFILT    Filter that excludes reads without the specified base(s) in the indent region
# --blist BARCODELIST   List of correct barcodes
#
# Sequence design:
# [Indent][Barcode][Spacer][DNA sequence]
# [AT][NN][GCT][DNA]
#
# Inputs:
# 1. A FASTQ file containing first-in-pair (Read 1) reads
# 2. A FASTQ file containing second-in-pair (Read 2) reads
#
# Outputs:
# 1. A Read 1 FASTQ file with barcodes added to the FASTQ header
# 2. A Read 2 FASTQ file with barcodes added to the FASTQ header
# 3. A text file summarizing barcode stats
#
###############################################################

################
#    Modules   #
################
from argparse import ArgumentParser
from itertools import zip_longest
import pandas as pd
import numpy as np
import sys


###############################
#        Main Function        #
###############################
def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--read1", action="store", dest="read1", type=str,
                        help="Input FASTQ file for Read 1 (unzipped)", required=True)
    parser.add_argument("--read2", action="store", dest="read2", type=str,
                        help="Input FASTQ file for Read 2 (unzipped)", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="Output SSCS BAM file", type=str,
                        required=True)
    parser.add_argument("--barcode", action="store", dest="barcode",
                        help="Barcode pattern (N = random/barcode, X = fixed/spacer, A|C|G|T = filtered bases) \n"
                             "e.g. XXNNXX means barcode flanked by 2 spacers on either side \n"
                             "e.g. NNGT means barcode is followed by two spacers matching 'GT'",
                        type=str, required=False)
    # parser.add_argument("--blen", action="store", dest="blen", help="Barcode length", type=int, required=True)
    # parser.add_argument("--slen", action="store", dest="slen", help="Spacer length (sequence between barcode and DNA)",
    #                     type=int, required=True)
    # parser.add_argument("--sfilt", action="store", dest="sfilt", type=str,
    #                     help="Spacer filter that excludes reads without the specified base(s) in the spacer region",
    #                     required=False)
    # parser.add_argument("--ilen", action="store", dest="ilen",
    #                     help="Indent length (sequence between adapter and barcode)", type=int, required=False)
    # parser.add_argument("--ifilt", action="store", dest="ifilt", type=str,
    #                     help="Indent filter that excludes reads without the specified base(s) in the indent region",
    #                     required=False)
    parser.add_argument("--blist", action="store", dest="blist", type=str, help="List of correct barcodes",
                        required=False)
    # parser.add_argument("--blistp", action="store", dest="blistp", type=str,
    #                     help="Barcode pattern for list of barcodes (e.g. NNXX means 2 barcode bases followed by 2 spacers",
    #                     required=False)
    args = parser.parse_args()

    ######################
    #       SETUP        #
    ######################
    # === Initialize input and output files ===
    read1 = open(args.read1, "r")
    read2 = open(args.read2, "r")
    r1_output = open('{}_barcode_R1.fastq'.format(args.outfile), "w")
    r2_output = open('{}_barcode_R2.fastq'.format(args.outfile), "w")
    stats = open('{}/barcode_stats.txt'.format(args.outfile.rsplit(sep="/", maxsplit=1)[0]), 'a')

    # === Initialize counters ===
    readpair_count = 0
    nospacer = 0
    noindent = 0
    bad_barcode = 0
    good_barcode = 0

    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    # Column in the following  corresponds to A, C, G, T, N
    r1_spacer_counter = pd.DataFrame(0, index=np.arange(args.slen), columns=nuc_lst)
    r2_spacer_counter = pd.DataFrame(0, index=np.arange(args.slen), columns=nuc_lst)
    r1_base_counter = pd.DataFrame(0, index=np.arange(args.blen), columns=nuc_lst)
    r2_base_counter = pd.DataFrame(0, index=np.arange(args.blen), columns=nuc_lst)

    # Rename rows
    r1_spacer_counter.index.names = ['R1_spacer']
    r2_spacer_counter.index.names = ['R2_spacer']
    r1_base_counter.index.names = ['R1_barcode']
    r2_base_counter.index.names = ['R2_barcode']

    # Check if there's indents
    if args.ilen is not None:
        r1_indent_counter = pd.DataFrame(0, index=np.arange(args.ilen), columns=nuc_lst)
        r2_indent_counter = pd.DataFrame(0, index=np.arange(args.ilen), columns=nuc_lst)
        r1_indent_counter.index_names = ['R1_indent']
        r2_indent_counter.index_names = ['R2_indent']

    ######################
    #  Extract barcodes  #
    ######################
    for r1, r2 in zip(zip_longest(*[read1] * 4), zip_longest(*[read2] * 4)):
        readpair_count += 1

        # Remove new line '\n' from str and separate using variables
        r1_header = r1[0].rstrip()
        r1_seq = r1[1].rstrip()
        r1_qual = r1[3].rstrip()

        r2_header = r2[0].rstrip()
        r2_seq = r2[1].rstrip()
        r2_qual = r2[3].rstrip()

        # Isolate indent
        if args.ilen is not None:
            r1_indent = r1_seq[0:args.ilen]
            r2_indent = r2_seq[0:args.ilen]

            r1_seq = r1_seq[args.ilen:len(r1_seq)]
            r2_seq = r2_seq[args.ilen:len(r2_seq)]

            r1_qual = r1_qual[args.ilen:len(r1_qual)]
            r2_qual = r2_qual[args.ilen:len(r2_qual)]

            # Count indent bases
            for i in range(len(r1_indent)):
                r1_indent_counter.iloc[i, nuc_lst.index(r1_indent[i])] += 1
                r2_indent_counter.iloc[i, nuc_lst.index(r2_indent[i])] += 1

        # Isolate spacer
        r1_spacer = r1_seq[args.blen:args.blen + args.slen]
        r2_spacer = r2_seq[args.blen:args.blen + args.slen]

        # Count spacer bases
        for i in range(len(r1_spacer)):
            r1_spacer_counter.iloc[i, nuc_lst.index(r1_spacer[i])] += 1
            r2_spacer_counter.iloc[i, nuc_lst.index(r2_spacer[i])] += 1

        # Check spacer filter
        if args.sfilt is not None and (r1_spacer != args.sfilt or r2_spacer != args.sfilt):
            nospacer += 1
        # Check indent filter
        elif args.ifilt is not None and (r1_indent != args.ifilt or r2_indent != args.ifilt):
            noindent += 1
        else:
            # Isolate barcodes
            r1_barcode = r1_seq[:args.blen]
            r2_barcode = r2_seq[:args.blen]

            # Count barcode bases
            for i in range(len(r1_barcode)):
                r1_base_counter.iloc[i, nuc_lst.index(r1_barcode[i])] += 1
                r2_base_counter.iloc[i, nuc_lst.index(r2_barcode[i])] += 1

            if r1_barcode.count("N") == 0 and r2_barcode.count("N") == 0:
                good_barcode += 1
                # Extract barcode from sequence and quality scores
                r1_seq = r1_seq[args.blen + args.slen:]
                r2_seq = r2_seq[args.blen + args.slen:]

                r1_qual = r1_qual[args.blen + args.slen:]
                r2_qual = r2_qual[args.blen + args.slen:]

                # Add barcode and read number to header
                r1_header = '{}|{}{}/{}'.format(r1_header.split(" ")[0], r1_barcode, r2_barcode, "1")
                r2_header = '{}|{}{}/{}'.format(r2_header.split(" ")[0], r1_barcode, r2_barcode, "2")

                # Write read to output file
                r1_output.write('{}\n{}\n+\n{}\n'.format(r1_header, r1_seq, r1_qual))
                r2_output.write('{}\n{}\n+\n{}\n'.format(r2_header, r2_seq, r2_qual))

            else:
                bad_barcode += 1

    r1_output.close()
    r2_output.close()

    # System output
    sys.stderr.write("Total sequences: {}\n".format(readpair_count))
    sys.stderr.write("Missing spacer: {}\n".format(nospacer))
    sys.stderr.write("Missing indent: {}\n".format(noindent))
    sys.stderr.write("Bad barcodes: {}\n".format(bad_barcode))
    sys.stderr.write("Passing barcodes: {}\n".format(good_barcode))

    # Output stats file
    stats.write("##########\n{}\n##########".format(args.outfile.split(sep="/")[-1]))

    if args.ilen is not None:
        stats.write(
            '\nTotal sequences: {}\nMissing indent: {}\nBad barcodes: {}\nPassing barcodes: {}\nMissing spacer: {}\n'.format(readpair_count,
                                                                                                         noindent,
                                                                                                         bad_barcode,
                                                                                                         good_barcode,
                                                                                                         nospacer))
        # Evaluate stats as %
        stats.write('---INDENT---\n{}\n-----------\n{}\n\n'.format(r1_indent_counter.apply(lambda x: x / x.sum(), axis=1),
                                                                             r2_indent_counter.apply(lambda x: x / x.sum(), axis=1)))

    else:
        stats.write('\nTotal sequences: {}\nMissing spacer: {}\nBad barcodes: {}\nPassing barcodes: {}\n'.format(readpair_count,
                                                                                                             nospacer,
                                                                                                             bad_barcode,
                                                                                                             good_barcode))

    stats.write('---BARCODE---\n{}\n-----------\n{}\n'.format(r1_base_counter.apply(lambda x: x / x.sum(), axis=1),
                                                              r2_base_counter.apply(lambda x: x / x.sum(), axis=1)))
    stats.write('---SPACER---\n{}\n-----------\n{}\n\n'.format(r1_spacer_counter.apply(lambda x: x / x.sum(), axis=1),
                                                               r2_spacer_counter.apply(lambda x: x / x.sum(), axis=1)))

    stats.close()


if __name__ == "__main__":
    main()






