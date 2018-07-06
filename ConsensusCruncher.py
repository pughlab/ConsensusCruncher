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
import argparse
# import extract_barcodes
import configparser
# from configobj import ConfigObj

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
    # Create directory for barcode extracted FASTQ files and BAM files
    fastq_dir = '{}/fastq_tag'.format(args.f_output)
    bam_dir = '{}/bamfiles'.format(args.f_output)
    # Check if dir exists and there's permission to write
    if not os.path.exists(fastq_dir) and not os.path.exists(bam_dir) and os.access(args.f_output, os.W_OK):
        os.makedirs(fastq_dir)
        os.makedirs(bam_dir)

    # Setup variables
    # if args.blist is not None and args.bpattern is not None:






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
    sub = main_p.add_subparsers(help='sub-command help', dest='subparser_name')
    sub.required = True

    #############
    # fastq2bam #
    #############
    # Help messages
    mode_fastq2bam_help = "Extract molecular barcodes from paired-end sequencing reads using a barcode list and/or " \
                          "a barcode pattern."
    # Set args
    sub_a = sub.add_parser('fastq2bam', help=mode_fastq2bam_help)
    sub_a.add_argument('-c', '--config', metavar="FILE", type=str, default=None,
                       help="Specify config file. Commandline option overrides config file (Use config template).")

    args, remaining_args = sub_a.parse_known_args()

    defaults = {"fastqs":"default"}

    if args.config:
        config = configparser.ConfigParser()
        config.read([args.config])
        defaults.update(dict(config.items("fastq2bam")))

        # sub_a = argparse.ArgumentParser(parents=[config])
        # sub_a.set_defaults(**defaults)

    # Parse commandline arguments
    sub_a.add_argument('-f', '--fastqs', dest='f_input', metavar="FASTQ", type=str, nargs=2,
                        help='Two paired-end FASTQ files.')
    sub_a.add_argument('-o', '--output', dest='f_output', metavar="OUTPUT_DIR", type=str,
                       help="Output directory, where barcode extracted FASTQ and BAM files will be placed in "
                            "subdirectories 'fastq_tag' and 'bamfiles' respectively (dir will be created if they "
                            "do not exist).")
    sub_a.add_argument('-r', '--ref', metavar="REF", help='Reference (BWA index).', type=str)
    sub_a.add_argument('-b', '--bpattern', metavar="BARCODE_PATTERN", type=str, default=None,
                        help="Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases).")
    sub_a.add_argument('-l', '--blist', metavar="BARCODE_LIST", type=str, default=None,
                       help='List of barcodes (Text file with unique barcodes on each line).')
    sub_a.set_defaults(func=fastq2bam)

    if args.config:
        sub_a.parse_args(remaining_args)

    #############
    # consensus #
    #############
    # Help messages
    bedfile_help = "Bedfile, default: cytoBand.txt. " \
                   "WARNING: It is HIGHLY RECOMMENDED that you use the default cytoBand.txt and" \
                   "not to include your own bedfile. This option is mainly intended for non-human" \
                   "genomes, where a separate bedfile is needed for data segmentation. If you do" \
                   "choose to use your own bedfile, please format with the bed_separator.R tool.\n" \
                   "For small or non-human genomes where cytobands cannot be used for segmenting the"\
                   "data set, you may choose to turn off this option with '-b OFF' and process the" \
                   "data all at once (Division of data is only required for large data sets to offload" \
                   "the memory burden)."
    consensus_help = "Almalgamate duplicate reads in BAM files into single-strand consensus sequences (SSCS) and " \
                     "duplex consensus sequences (DCS). Single reads with complementary duplex strands can also be " \
                     "corrected with 'Singleton Correction'."

    # Determine code directory
    code_dir = os.path.join(os.path.dirname(__file__))
    print(code_dir)

    # Set args
    sub_b = sub.add_parser('consensus', help=consensus_help)
    sub_b.add_argument('-i', '--input', dest='c_input', help='Input directory.', required=True, type=str)
    sub_b.add_argument('-o', '--output', dest='c_output', required=True, type=str,
                       help="Output project directory for new files and folders to be created.")
    sub_b.add_argument('-s', '--scorrect', help="Singleton correction, default: True.", default=True,
                       choices=[True, False], type=bool)
    sub_b.add_argument('-b', '--bedfile', help=bedfile_help, default='{}/cytoband.txt'.format(code_dir), type=str)
    sub_b.add_argument('-c', '--cutoff', default=0.7, type=float,
                       help="Consensus cut-off, default: 0.7 (70%% of reads must have the same base to form a "
                            "consensus).")
    sub_b.add_argument('--cleanup', help="Remove intermediate files.")
    sub_b.set_defaults(func=consensus)

    # Parse args
    args = main_p.parse_args()

    if args.config is None and (args.f_input is None or args.f_output is None or args.ref is None):
        sub_a.error("Command line arguments must be provided if config file is not present.")

    # Check if either barcode pattern or list is set. At least one must be provided.
    if args.bpattern is None and args.blist is None:
        sub_a.error("At least one of -b or -l required.")
    # Check proper barcode design provided for barcode pattern
    elif re.findall(r'[^A|C|G|T|N]', args.bpattern):
        raise ValueError("Invalid barcode pattern containing characters other than A, C, G, T, and N.")
    # Check list for faulty barcodes in list
    elif args.blist is not None:
        blist = open(args.blist, "r").read().splitlines()
        if re.search("[^ACGTN]", "".join(blist)) is not None:
            raise ValueError("List contains invalid barcodes. Please specify barcodes with A|C|G|T.")
        else:
            args.func(args)
    else:
        args.func(args)



