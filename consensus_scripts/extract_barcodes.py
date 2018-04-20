#!/usr/bin/env python3

###############################################################
#
#                       Barcode Extractor
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
# --read1 READ1       Input FASTQ file for Read 1 (unzipped)
# --read2 READ2       Input FASTQ file for Read 2 (unzipped)
# --outfile OUTFILE   Output FASTQ files for Read 1 and Read 2 using given filename
# --blen BARCODELEN   Barcode length
# --slen SPACERLEN    Spacer length (This region is removed and only the barcode will be added to the header)
# --sfilt SPACERFILT  Filter that excludes reads without the specified base(s) in the spacer region
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
    parser.add_argument("--blen", action="store", dest="blen", help="Barcode length", type=int, required=True)
    parser.add_argument("--slen", action="store", dest="slen", help="Spacer length", type=int, required=True)
    parser.add_argument("--sfilt", action="store", dest="sfilt", type=str,
                        help="Spacer filter that excludes reads without the specified base(s) in the spacer region",
                        required=False)
    args = parser.parse_args()

    ######################
    #       SETUP        #
    ######################
    # === Initialize input and output files ===
    read1 = open(args.read1, "r")
    read2 = open(args.read2, "r")
    r1_output = open('{}_barcode_R1.fastq'.format(args.outfile), "w")
    r2_output = open('{}_barcode_R2.fastq'.format(args.outfile), "w")
    stats = open('{}.stats.txt'.format(args.outfile), 'a')

    # === Initialize counters ===
    readpair_count = 0
    nospacer = 0
    bad_barcode = 0
    good_barcode = 0

    nuc_lst = ['A', 'C', 'G', 'T', 'N']
    # Positions in the following lists corresponds to A, C, G, T, N
    r1_base_counter = [0, 0, 0, 0, 0]
    r1_spacer_counter = [0, 0, 0, 0, 0]
    r2_base_counter = [0, 0, 0, 0, 0]
    r2_spacer_counter = [0, 0, 0, 0, 0]

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

        # Isolate spacer
        r1_spacer = r1_seq[args.blen:args.blen + args.slen]
        r2_spacer = r2_seq[args.blen:args.blen + args.slen]

        r1_spacer_counter[nuc_lst.index(r1_spacer)] += 1
        r2_spacer_counter[nuc_lst.index(r2_spacer)] += 1

        # Check spacer filter
        if args.sfilt is not None and r1_spacer is not args.sfilt and r2_spacer is not args.sfilt:
            nospacer += 1
        else:
            # Isolate barcodes
            r1_barcode = r1_seq[:args.blen]
            r2_barcode = r2_seq[:args.blen]

            for base in list(r1_barcode):
                r1_base_counter[nuc_lst.index(base)] += 1

            for base in list(r2_barcode):
                r2_base_counter[nuc_lst.index(base)] += 1

            if r1_barcode.count("N") == 0 and r2_barcode.count("N") == 0:
                good_barcode += 1
                # Extract barcode from sequence and quality scores
                r1_seq = r1_seq[args.blen + args.slen :]
                r2_seq = r2_seq[args.blen + args.slen :]

                r1_qual = r1_qual[args.blen + args.slen :]
                r2_qual = r2_qual[args.blen + args.slen :]

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
    sys.stderr.write("Bad barcodes: {}\n".format(bad_barcode))
    sys.stderr.write("Passing barcodes: {}\n".format(good_barcode))

    # Output stats file
    stats.write(args.outfile)
    stats.write('\nTotal sequences: {}\nMissing spacer: {}\nBad barcodes: {}\nPassing barcodes: {}\n'.format(readpair_count,
                                                                                                         nospacer,
                                                                                                         bad_barcode,
                                                                                                         good_barcode))
    stats.write('R1 barcode (A, C, G, T, N): {}\nR1 spacer (A, C, G, T, N): {}\n'.format(r1_base_counter,
                                                                                       r1_spacer_counter))
    stats.write('R1 barcode (A, C, G, T, N): {}\nR1 spacer (A, C, G, T, N): {}\n\n'.format(r2_base_counter,
                                                                                       r2_spacer_counter))

    stats.close()


if __name__ == "__main__":
    main()






