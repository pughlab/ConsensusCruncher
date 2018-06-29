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
#  WARNING: Please change qsub parameters according to your cluster commands; default
#  qsub -q highmem.q SCRIPT
#
# ===================================================================================

import os
import sys


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
    subp_p.add_argument('-i', dest='input', help='Input directory.', required=True, type=str)
    subp_p.add_argument('-o', dest='output', help='Output project directory for new folders and files.', required=True,
                        type=str)
    subp_p.add_argument('-b', dest='bpattern', type=str,
                        help="Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases).")
    subp_p.add_argument('-l', dest='blist', type=str,
                        help='List of barcodes (Text file with unique barcodes on each line).')
    subp_p.set_defaults(func=fastq2bam)

    #############
    # consensus #
    #############
    # Help messages

    # Set args
    subp_ = subp.add_parser('consensus', help=)


    # Parse args
    args = main_p.parse_args()
    args.func(args)

