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
import matplotlib.pyplot as plt


#######################
#    Helper Function    #
#######################
def check_overlap(blist):
    """(list) -> bool
    Return boolean indicating whether or not there's overlapping barcodes within the list.

    >>> check_overlap(['AACT', 'AGCT'])
    False
    >>> check_overlap(['AACTCT', 'AACT'])
    True
    """
    blist = list(set(blist))
    print(blist)
    overlap = False
    for barcode in blist:
        print("Barcode is", barcode)
        over = [i for i in blist if barcode in i]
        
        if len(over) == 1:
            print("No other barcode contains ", barcode, ". No Overlap")
        if len(over) > 1:
            print(barcode, "Is contained in ", over)
            
            m= (min(over, key=len))
            print(m)
            over.remove(m)
            print(over)
            for item in over:
                if item[:len(m)] == m:
                    overlap= True
                    print( "There is overlapping barcodes")
                else:
                    print(item, "Does not start with ", barcode, ". No Overlap")
 


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
    # Raise error if neither a barcode list or pattern is provided
    if args.blist is None and args.bpattern is None:
        raise ValueError("No barcode specifications inputted. Please specify barcode list or pattern.")
    # == Barcode Pattern ==
    elif args.bpattern is not None:
        # Ensure valid barcode pattern provided
        if re.search("[^ACGTN]", args.bpattern) is not None:
            raise ValueError("Invalid barcode pattern inputted. Please specify pattern with A|C|G|T = fixed, "
                             "N = variable (e.g. 'ATNNGCT').")
        else:
            plen = len(args.bpattern)  # Pattern length
            b_index = list(find_all(args.bpattern, 'N'))  # Index of random barcode bases
            s_index = [x for x in list(range(0, plen)) if x not in b_index]  # Index of constant spacer bases
            spacer = ''.join([args.bpattern[x] for x in s_index])

            # Column in the following corresponds to A, C, G, T, N
            nuc_dict = create_nuc_dict(nuc_lst)
            r1_barcode_counter = pd.DataFrame(0, index=np.arange(plen), columns=nuc_lst)
            r2_barcode_counter = pd.DataFrame(0, index=np.arange(plen), columns=nuc_lst)
    # == Barcode list ==
    else:
        blist = open(args.blist, "r").read().splitlines()

        # Check list for faulty barcodes
        if re.search("[^ACGTN]", "".join(blist)) is not None:
            raise ValueError("Invalid barcode list inputted. Please specify barcodes with A|C|G|T.")
        # Check barcodes end with spacer T (necessary for T-tailed adapters and 3' dA overhang on the fragmented DNA sample)
        elif [s for s in blist if not s.endswith("T")] != []:
            raise ValueError("There is one or more barcodes in the list that do not end with 'T'.")
        # Check barcodes in list are not overlapping as its indicative of faulty design
        # (Difficult to differentiate whether the shorter or longer barcode is correct)
        elif check_overlap(blist):
            raise ValueError("There are overlapping barcodes in the list (difficult to determine which barcode is correct).")
        else:
            # Barcode counter: create dictionary with barcodes as keys and values as 0
            # - Barcodes may be of different lengths, so a tally of each barcode occurrence
            # is better than moderating the frequency of nuc bases at each barcode position
            r1_tag_dict = dict.fromkeys(blist, 0)
            r2_tag_dict = dict.fromkeys(blist, 0)

            # Write bad_barcodes to file
            r1_bad_barcodes = open('{}_r1_bad_barcodes.txt'.format(args.outfile), 'w')
            r2_bad_barcodes = open('{}_r2_bad_barcodes.txt'.format(args.outfile), 'w')


    ######################
    #  Extract barcodes  #
    ######################
    for r1, r2 in zip(read1, read2):
        readpair_count += 1

        # Check if R1 and R2 matches
        assert r1.id == r2.id

        # === Barcode pattern ===
        if args.bpattern is not None:
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
                r1_bc = ''.join([r1_barcode[x] for x in b_index])
                r2_bc = ''.join([r2_barcode[x] for x in b_index])

                r1.id = '{}|{}.{}/{}'.format(r1.id.split(" ")[0], r1_bc, r2_bc, "1")
                r2.id = '{}|{}.{}/{}'.format(r2.id.split(" ")[0], r1_bc, r2_bc, "2")
                r1.description = r1.id
                r2.description = r2.id

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

        # === Barcode list ===
        else:
            # Identify unique lengths of barcodes
            barcode_len = list(set(len(n) for n in blist))
            # Sort length from highest to lowest in case shorter barcodes overlap longer ones
            barcode_len.sort(reverse=True)

            # Check to see if both R1 and R2 barcodes are present in blist
            r1_status = False
            r2_status = False

            # Iterate through barcodes of different lengths
            for blen in barcode_len:
                # Remove new line '\n' from str and separate using barcode length
                r1_read, r1_barcode = extract_barcode(r1, blen)
                r2_read, r2_barcode = extract_barcode(r2, blen)

                # Check to see if barcode is valid
                if re.search("[^ACGT]", r1_barcode) is not None or re.search("[^ACGT]", r2_barcode) is not None:
                    bad_barcode += 1
                    if re.search("[^ACGT]", r1_barcode) is not None:
                        r1_bad_barcodes.write(r1_barcode + '\n')
                    if re.search("[^ACGT]", r2_barcode) is not None:
                        r2_bad_barcodes.write(r2_barcode + '\n')
                else:
                    # Save barcode and read if barcode is found in blist
                    if r1_barcode in blist:
                        r1_status = True
                        r1_r = r1_read
                        r1_b = r1_barcode[:blen-1]	# remove T from end of barcode

                    if r2_barcode in blist:
                        r2_status = True
                        r2_r = r2_read
                        r2_b = r2_barcode[:blen-1]  # remove T from end of barcode

            # If R1 and R2 barcodes are both valid
            if r1_status and r2_status:
                good_barcode += 1  # Note: number of barcodes is per paired reads

                # Add to barcode counter
                r1_tag_dict[r1_b+'T'] += 1
                r2_tag_dict[r2_b+'T'] += 1

                # Add barcode and read number to header of fastq
                r1_r.id = '{}|{}.{}/{}'.format(r1_r.id.split(" ")[0], r1_b, r2_b, "1")
                r2_r.id = '{}|{}.{}/{}'.format(r2_r.id.split(" ")[0], r1_b, r2_b, "2")
                # Update description so ID is not repeated twice in FASTQ header
                r1_r.description = r1_r.id
                r2_r.description = r2_r.id

                SeqIO.write(r1_r, r1_output, "fastq")
                SeqIO.write(r2_r, r2_output, "fastq")
            else:
                # Note bad barcodes always correspond to length of shortest barcode as we can't determine the original
                # barcode length
                bad_barcode += 1
                if r1_status == False:
                    r1_bad_barcodes.write(r1_barcode + '\n')
                if r2_status == False:
                    r2_bad_barcodes.write(r2_barcode + '\n')

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
    # == Barcode pattern ==
    if args.bpattern is not None:
        stats.write('---BARCODE---\n{}\n-----------\n{}\n'.format(r1_barcode_counter.apply(lambda x: x / x.sum(), axis=1),
                                                                  r2_barcode_counter.apply(lambda x: x / x.sum(), axis=1)))
    # == Barcode list ==
    else:
        # Convert dict to dataframes
        r1_df = pd.DataFrame(sorted(r1_tag_dict.items(), key=lambda kv: (len(kv[0]), kv[0])), columns=["Barcode", "R1_Count"])
        r2_df = pd.DataFrame(sorted(r2_tag_dict.items(), key=lambda kv: (len(kv[0]), kv[0])), columns=["Barcode", "R2_Count"])

        # Merge dataframes and sum total count
        df_merge = pd.merge(r1_df, r2_df, on="Barcode")
        df_merge['Total'] = df_merge['R1_Count'] + df_merge['R2_Count']

        # Order dataframe
        df_merge = df_merge.sort_values(by="Total", ascending=False)

        # Write stats to file
        stats.write('---BARCODE---\n{}\n'.format(df_merge))

        # == Create histogram for barcode stats ==
        fig, ax = plt.subplots()
        ax.set_xlim(0, len(df_merge.index))  # Set x-axis range to number of barcodes

        ind = np.arange(len(df_merge.index))  # the x locations (number of barcodes)
        width = 0.35  # the width of the bars
        p1 = ax.bar(ind, df_merge['R1_Count'], width, color='g')
        p2 = ax.bar(ind + width, df_merge['R2_Count'], width, color='y')

        # Label axis
        ax.set_xticks(ind + width)
        ax.set_xticklabels(df_merge['Barcode'])
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        plt.gcf().subplots_adjust(bottom=0.15)  # Add space to make sure x labels aren't cut off
        plt.tick_params(axis='x', which='both', top=False)
        plt.tick_params(axis='y', which='both', right=False)

        # Set legends and labels
        ax.legend((p1[0], p2[0]), ('Read1', 'Read2'))
        ax.set_title('Barcode frequency')
        plt.ylabel('Count')

        plt.savefig('{}_barcode_stats.png'.format(args.outfile))

    stats.close()

if __name__ == "__main__":
    main()
