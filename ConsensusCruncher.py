#!/usr/bin/env python3

# ===================================================================================
#
#  FILE:         ConsensusCruncher.py
#
#  USAGE:
#  ConsensusCruncher.py fastq2bam -i input_dir -o output_dir
#
#
#  ConsensusCruncher.py consensus -i input_dir -o output_dir
#
#
#  OPTIONS:
#
#    -i  Input bamfile directory [MANDATORY]
#    -o  Output project directory [MANDATORY]
#    -s  Singleton correction, default: ON (use "OFF" to disable)
#    -b  Bedfile, default: cytoBand.txt
#        WARNING: It is HIGHLY RECOMMENDED that you use the default cytoBand.txt and
#        not to include your own bedfile. This option is mainly intended for non-human
#        genomes, where a separate bedfile is needed for data segmentation. If you do
#        choose to use your own bedfile, please format with the bed_separator.R tool.
#
#        For small or non-human genomes where cytobands cannot be used for segmenting the
#        data set, you may choose to turn off this option with "-b OFF" and process the
#        data all at once (Division of data is only required for large data sets to offload
#        the memory burden).
#
#    -c  Consensus cut-off, default: 0.7 (70% of reads must have the same base to form
#        a consensus)
#    -q  qusb directory, default: output/qsub
#    -h  Show this message
#
#  DESCRIPTION:
#
#  This script amalgamates duplicate reads in bamfiles into single-strand consensus
#  sequences (SSCS), which are subsequently combined into duplex consensus sequences
#  (DCS). Singletons (reads lacking duplicate sequences) are corrected, combined
#  with SSCS to form SSCS + SC, and further collapsed to form DCS + SC. Finally,
#  files containing all unique molecules (a.k.a. no duplicates) are created for SSCS
#  and DCS.
#
#  Note: Script will create a "consensus" directory under the project directory and
#  sub-directories corresponding to each bamfile in the input directory.
#
# ===================================================================================

import os
import sys
import re

####################
# Helper functions #
####################
def check_barcode(bpattern):
    """
    Check barcodes only contain A, C, G, T, and N bases.

    :param bpattern(str): the sequence of a barcode, e.g. 'ATNNCGT'
    :return:
    """


def fastq2bam(args):
    """
    Extract molecular barcodes from paired-end sequencing reads using a barcode list,
    pattern, or the two combined. Remove constant spacer bases and combine paired
    barcodes before adding to the header of each read in FASTQ files.

    Barcode-extracted FASTQ files are written to the 'fastq_tag' directory and are
    subsequenntly aligned with BWA. BAM files are written to a 'bamfile' directory
    under the specified project folder.

    BARCODE DESIGN:
    You can input either a barcode list or barcode pattern or both. If both are provided, barcodes will first be matched
    with the list and then the constant spacer bases will be removed before the barcode is added to the header.

    N = random / barcode bases
    A | C | G | T = constant spacer bases
    e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front and 'GT' behind.


    :param args:

    """
    # Check if either barcode pattern or list is set. At least one must be provided.


    # Check proper barcode design provided for barcode pattern
    try:
        if re.findall(r'[^A|C|G|T|N]', args.bpattern):
            raise ValueError
    except ValueError:
        print("Invalid barcode pattern containing characters other than A, C, G, T, and N.")

    # Create directory (check if dir exists and permissions to write)
    fastq_dir = '{}/fastq_tag'.format(args.f_output)
    bam_dir = '{}/bamfiles'.format(args.f_output)
    if not os.path.exists(fastq_dir) and not os.path.exists(bam_dir) and os.access(args.f_output, os.W_OK):
        os.makedirs(fastq_dir)
        os.makedirs(bam_dir)





def consensus(args):
    """
    Using unique molecular identifiers (UMIs), duplicate reads from the same molecule are amalgamated into single-strand
    consensus sequences (SSCS). If complementary strands are present, SSCSs can be subsequently combined to form duplex
    consensus sequences (DCS).

    If 'Singleton Correction' (SC) is enabled, single reads (singletons) can be error suppressed using complementary
    strand. These corrected singletons can be merged with SSCSs to be further collapsed into DCSs + SC.

    Finally, a BAM file containing only unique molecules (i.e. no duplicates) is created by merging DCSs, remaining
    SSCSs (those that could not form DCSs), and remaining singletons (those that could not be corrected).


    :param args:

    """



if __name__ == '__main__':
    # Mode parser
    main_p = argparse.ArgumentParser()
    subp = main_p.add_subparsers(help='sub-command help')

    #############
    # fastq2bam #
    #############
    # Help messages
    mode_fastq2bam_help = "Extract molecular barcodes from paired-end sequencing reads using a barcode list" \
                          "pattern, or the two combined. Remove constant spacer bases and combine paired" \
                          "barcodes before adding to the header of each read in FASTQ files.\n" \
                          "Barcode-extracted FASTQ files are written to the 'fastq_tag' directory and are" \
                          "subsequenntly aligned with BWA. BAM files are written to a 'bamfile' directory" \
                          "under the specified project folder.\n" \
                          "BARCODE DESIGN:" \
                          "You can input either a barcode list or barcode pattern or both. If both are provided," \
                          "barcodes will first be matched with the list and then the constant spacer bases will" \
                          "be removed before the barcode is added to the header.\n" \
                          "N = random barcode bases" \
                          " A | C | G | T = constant spacer bases" \
                          "e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front and 'GT' behind."

    # Set args
    subp_p = subp.add_parser('fastq2bam', help=mode_fastq2bam_help)
    subp_p.add_argument('-i', dest='f_input', required=True, type=str,
                        help='Input FASTQ directory (Note: filenames must contain R1/R2 to differentiate paired-end '
                             'reads).')
    subp_p.add_argument('-o', dest='f_output', help='Output project directory for new folders and files to be created.',
                        required=True, type=str)
    subp_p.add_argument('-r', dest='ref', help='Reference (BWA index).', required=TRUE, type=str)
    subp_p.add_argument('-b', dest='bpattern', type=str,
                        help="Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases).")
    subp_p.add_argument('-l', dest='blist', type=str,
                        help='List of barcodes (Text file with unique barcodes on each line).')
    subp_p.set_defaults(func=fastq2bam)


    #############
    # consensus #
    #############
    # Help messages
    bedfile_help = "Bedfile, default: cytoBand.txt" \
                   "WARNING: It is HIGHLY RECOMMENDED that you use the default cytoBand.txt and" \
                   "not to include your own bedfile. This option is mainly intended for non-human" \
                   "genomes, where a separate bedfile is needed for data segmentation. If you do" \
                   "choose to use your own bedfile, please format with the bed_separator.R tool.\n" \
                   "For small or non-human genomes where cytobands cannot be used for segmenting the"\
                   "data set, you may choose to turn off this option with '-b OFF' and process the" \
                   "data all at once (Division of data is only required for large data sets to offload" \
                   "the memory burden)."
    # Determine code directory
    code_dir

    # Set args
    subp_p = subp.add_parser('consensus', help=)
    subp_p.add_argument('-i', dest='c_input', help='Input directory.', required=True, type = str)
    subp_p.add_argument('-o', dest='c_output', help="Output project directory for new files and folders to be created.",
                        required=True, type=str)
    subp_p.add_argument('-s', dest='singcor', help="Singleton correction, default: True.",
                        default=True, action='store', choices=[True, False], type=bool)
    subp_p.add_argument('-b', dest='bedfile', help=bedfile_help, default='{}/cytoband.txt'.format(code_dir), type=str)
    subp_p.add_argument('-c', dest='cutoff', default=0.7, type=float,
                        help="Consensus cut-off, default: 0.7 (70% of reads must have the same base to form a "
                             "consensus).")
    subp_p.set_defaults(func=consensus)

    # Parse args
    args = main_p.parse_args()

    if args.bpattern is None and args.blist is None:
        args.error("At least one of -b or -l required.")
    else:
        args.func(args)

