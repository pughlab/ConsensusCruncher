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
# --bpattern BPATTERN   Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases)
# --blist BARCODELIST   List of correct barcodes
#
# Barcode design:
# N = random / barcode bases
# A | C | G | T = Fixed spacer bases
# e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front and 'GT' behind
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
from gzip import open as gzopen
from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import sys


#######################
#    Helper Function    #
#######################
def find_all(a_str, sub):
    """(str, str) -> int
    Return index of substring in string.
    """
    sub_index=[i for i, x in enumerate(list(a_str)) if x==sub]    
    return sub_index


def create_nuc_dict(nuc_lst):
    """ (list) -> dict
    Takes the nucleotide list and converts it to a dictionary of binary arrays.
    e.g. A => {"A":[1, 0, 0, 0, 0]}
    """
    nuc_dict = {}
    for nuc_x in nuc_lst:
        nuc_arr = np.array([nuc_y == nuc_x for nuc_y in nuc_lst]).astype(int)
        nuc_dict.update({nuc_x:nuc_arr})
    return nuc_dict


def seq_to_mat(seq, nuc_dict):
    """ (string, dict) -> np.matrix
    Converts a string of nucleotides to a matrix of binary nucleotides.
    """
    nuc_map = [nuc_dict[each_nuc] for each_nuc in list(seq)]
    nuc_mat = np.asmatrix(nuc_map)
    return nuc_mat


def extract_barcode(read, plen):
    """
    Extract barcode from Seq and Phred quality.

    :param read: A SeqIO object.
    :type read: object
    :param plen: The length of the barcode.
    :type plen: num
    :returns: A SeqIO object with barcode removed and a barcode string.
    """
    barcode = str(read.seq[:plen])
    read_b = read[plen:]

    return read_b, barcode


#######################
#    Main Function    #
#######################
def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--read1", action="store", dest="read1", type=str,
                        help="Input FASTQ file for Read 1 (unzipped)", required=True)
    parser.add_argument("--read2", action="store", dest="read2", type=str,
                        help="Input FASTQ file for Read 2 (unzipped)", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="Absolute path to output SSCS BAM file",
                        type=str, required=True)
    parser.add_argument("--bpattern", action="store", dest="bpattern", type=str, required=False, default=None,
                        help="Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases) \n"
                             "e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front, "
                             "followed by 'GT' \n")
    parser.add_argument("--blist", action="store", dest="blist", type=str, help="List of correct barcodes",
                        default=None, required=False)
    args = parser.parse_args()

    ######################
    #       SETUP        #
    ######################
    # === Initialize input and output files ===
    # Check if file is zipped
    if 'gz' in args.read1:
        read1 = SeqIO.parse(gzopen(args.read1, "rt"), "fastq")
        read2 = SeqIO.parse(gzopen(args.read2, "rt"), "fastq")
    else:
        read1 = SeqIO.parse(open(args.read1, "rU"), "fastq")
        read2 = SeqIO.parse(open(args.read2, "rU"), "fastq")

    r1_output = open('{}_barcode_R1.fastq'.format(args.outfile), "w")
    r2_output = open('{}_barcode_R2.fastq'.format(args.outfile), "w")
    stats = open('{}_barcode_stats.txt'.format(args.outfile.rsplit(sep="/", maxsplit=1)[0]), 'a')
    
    # === Initialize counters ===
    readpair_count = 0
    bad_spacer = 0
    bad_barcode = 0
    good_barcode = 0
    
    nuc_lst = ['A', 'C', 'G', 'T', 'N']

    # === Define barcodes ===
    # Check if list of barcodes is provided
    if args.blist is not None:
        blist = open(args.blist, "r").read().splitlines()
        plen = len(blist[0])  # Length of barcode

        # Check list for faulty barcodes
        if re.search("[^ACGTN]", "".join(blist)) is not None:
            raise ValueError("Invalid barcode list inputted. Please specify barcodes with A|C|G|T.")
    
    # Check if barcode pattern is provided
    if args.bpattern is not None:
        # Ensure valid barcode pattern provided
        if re.search("[^ACGTN]", args.bpattern) is not None:
            raise ValueError("Invalid barcode pattern inputted. Please specify pattern with A|C|G|T = fixed, "
                             "N = variable (e.g. 'ATNNGCT').")

        plen = len(args.bpattern)  # Pattern length
        b_index = list(find_all(args.bpattern, 'N'))  # Index of random barcode bases
        s_index = [x for x in list(range(0, plen)) if x not in b_index]  # Index of constant spacer bases
        spacer = ''.join([args.bpattern[x] for x in s_index])

    # Raise error if neither a barcode list or pattern is provided
    if args.blist is None and args.bpattern is None:
        raise ValueError("No barcode specifications inputted. Please specify barcode list or pattern.")

    # Column in the following corresponds to A, C, G, T, N
    nuc_dict = create_nuc_dict(nuc_lst)
    r1_barcode_counter = pd.DataFrame(0, index=np.arange(plen), columns=nuc_lst)
    r2_barcode_counter = pd.DataFrame(0, index=np.arange(plen), columns=nuc_lst)

    ######################
    #  Extract barcodes  #
    ######################
    for r1, r2 in zip(read1, read2):
        readpair_count += 1

        # Check if R1 and R2 matches
        assert r1.id == r2.id

        # Remove new line '\n' from str and separate using variables
        r1, r1_barcode = extract_barcode(r1, plen)
        r2, r2_barcode = extract_barcode(r2, plen)

        # Check to see if barcode is valid
        if re.search("[^ACGT]", r1_barcode) is not None or re.search("[^ACGT]", r2_barcode) is not None:
            bad_barcode += 1
        else:
            # Count barcode bases
            r1_barcode_counter += seq_to_mat(r1_barcode, nuc_dict)
            r2_barcode_counter += seq_to_mat(r2_barcode, nuc_dict)

            # Add barcode and read number to header
            if args.bpattern is not None:
                r1_bc = ''.join([r1_barcode[x] for x in b_index])
                r2_bc = ''.join([r2_barcode[x] for x in b_index])

            r1.id = '{}|{}{}/{}'.format(r1.id.split(" ")[0], r1_bc, r2_bc, "1")
            r2.id = '{}|{}{}/{}'.format(r2.id.split(" ")[0], r1_bc, r2_bc, "2")
            r1.description = r1.id
            r2.description = r2.id

            # Isolate barcode from sequence
            if args.blist is not None:
                if r1_barcode in blist and r2_barcode in blist:
                    good_barcode += 1
                    # Write read to output file
                    SeqIO.write(r1, r1_output, "fastq")
                    SeqIO.write(r2, r2_output, "fastq")
                else:
                    bad_barcode += 1

            else:
                r1_spacer = ''.join([r1_barcode[x] for x in s_index])
                r2_spacer = ''.join([r2_barcode[x] for x in s_index])

                # Check if spacer is correct
                if r1_spacer == spacer and r2_spacer == spacer:
                    good_barcode += 1
                    # Write read to output file
                    SeqIO.write(r1, r1_output, "fastq")
                    SeqIO.write(r2, r2_output, "fastq")
                else:
                    bad_spacer += 1

    r1_output.close()
    r2_output.close()

    # System output
    sys.stderr.write("Total sequences: {}\n".format(readpair_count))
    sys.stderr.write("Missing spacer: {}\n".format(bad_spacer))
    sys.stderr.write("Bad barcodes: {}\n".format(bad_barcode))
    sys.stderr.write("Passing barcodes: {}\n".format(good_barcode))

    # Output stats file
    stats.write("##########\n{}\n##########".format(args.outfile.split(sep="/")[-1]))
    stats.write(
        '\nTotal sequences: {}\nMissing spacer: {}\nBad barcodes: {}\nPassing barcodes: {}\n'.format(readpair_count,
                                                                                                     bad_spacer,
                                                                                                     bad_barcode,
                                                                                                     good_barcode))
    stats.write('---BARCODE---\n{}\n-----------\n{}\n'.format(r1_barcode_counter.apply(lambda x: x / x.sum(), axis=1),
                                                              r2_barcode_counter.apply(lambda x: x / x.sum(), axis=1)))

    stats.close()


if __name__ == "__main__":
    main()
