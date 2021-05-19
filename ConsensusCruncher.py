#!/usr/bin/env python3

import os
import sys
import re
import argparse
import configparser
from subprocess import Popen, PIPE, call

def sort_index(bam, samtools):
    """
    Sort and index BAM file.

    :param bam: Path to BAM file.
    :type bam: str
    :param samtools: Path to samtools.
    :type samtools: str
    :returns: Path to sorted BAM file.
    """
    identifier = bam.split('.bam', 1)[0]
    sorted_bam = '{}.sorted.bam'.format(identifier)

    sam1 = Popen((samtools + ' view -bu ' + bam).split(' '), stdout=PIPE)
    sam2 = Popen((samtools + ' sort -').split(' '), stdin=sam1.stdout, stdout=open(sorted_bam, 'w'))
    sam2.communicate()
    os.remove(bam)
    call("{} index {}".format(samtools, sorted_bam).split(' '))

    return sorted_bam


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

    N = random / barcode base
    A | C | G | T = constant spacer bases
    e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front and 'GT' behind.
    """
    # Create directory for barcode extracted FASTQ files and BAM files
    fastq_dir = '{}/fastq_tag'.format(args.output)
    bam_dir = '{}/bamfiles'.format(args.output)

    # Check if dir exists and there's permission to write
    if not os.path.exists(fastq_dir) and os.access(args.output, os.W_OK):
        os.makedirs(fastq_dir)
    if not os.path.exists(bam_dir) and os.access(args.output, os.W_OK):
        os.makedirs(bam_dir)

    # Set file variables
    filename = os.path.basename(args.fastq1).split(args.name, 1)[0]
    outfile = "{}/{}".format(fastq_dir, filename)

    ####################
    # Extract barcodes #
    ####################
    if args.blist is not None and args.bpattern is not None:
        extractb_cmd = "{}/ConsensusCruncher/extract_barcodes.py --read1 {} --read2 {} --outfile {} --bpattern {} "
        "--blist {}".format(code_dir, args.fastq1, args.fastq2, outfile, args.bpattern, args.blist)
    elif args.blist is None:
        extractb_cmd = "{}/ConsensusCruncher/extract_barcodes.py --read1 {} --read2 {} --outfile {} --bpattern {}".format(
            code_dir, args.fastq1, args.fastq2, outfile, args.bpattern)
    else:
        extractb_cmd = "{}/ConsensusCruncher/extract_barcodes.py --read1 {} --read2 {} --outfile {} --blist {}".format(
            code_dir, args.fastq1, args.fastq2, outfile, args.blist)
    
    print(extractb_cmd)
    os.system(extractb_cmd)


    #Create directories for bad barcodes and barcode distribution histograms
    if args.blist is not None:
       bad_barcode_dir = '{}/fastq_tag/bad_barcode'.format(args.output)
       barcode_dist_dir = '{}/fastq_tag/barcode_dist'.format(args.output)

       if not os.path.exists(bad_barcode_dir) and os.access(args.output, os.W_OK):
           os.makedirs(bad_barcode_dir)

       if not os.path.exists(barcode_dist_dir) and os.access(args.output, os.W_OK):
           os.makedirs(barcode_dist_dir)

   #     Move files 
       os.rename('{}/{}_r1_bad_barcodes.txt'.format(fastq_dir, filename),
             '{}/{}_r1_bad_barcodes.txt'.format(bad_barcode_dir, filename))
       os.rename('{}/{}_r2_bad_barcodes.txt'.format(fastq_dir, filename),
             '{}/{}_r2_bad_barcodes.txt'.format(bad_barcode_dir, filename))
       os.rename('{}/{}_barcode_stats.png'.format(fastq_dir, filename),
             '{}/{}_barcode_stats.png'.format(barcode_dist_dir, filename))

    #############
    # BWA Align #
    #############
    # Command split into chunks and bwa_id retained as str repr
    picard =  'java -jar ' + args.picard + ' AddOrReplaceReadGroups' # "java -jar /mnt/work1/software/picard/2.10.9/picard.jar AddOrReplaceReadGroups"
    
    bwa_cmd = args.bwa + 'mem -M -t4'
    
    #bwa_id = "@RG\tID:1\tSM:" + filename + "\tPL:Illumina"
    bwa_args = '{} {}_barcode_R1.fastq {}_barcode_R2.fastq'.format(args.ref, outfile, outfile)
    
    bwa_cmd = args.bwa + ' mem -M -t4 ' + bwa_args
    print(bwa_cmd)
    bwa = Popen(bwa_cmd.split(' '), stdout=PIPE)
    #print(bwa)
    # # Sort BAM (BWA output piped into samtools for sorting before writing into bam)
    sam1 = Popen((args.samtools + ' view -bhS -').split(' '), stdin=bwa.stdout, stdout=PIPE)
    sam2 = Popen((args.samtools + ' sort -').split(' '), stdin=sam1.stdout,
                  stdout=open('{}/{}.sort.bam'.format(bam_dir, filename), 'w'))
    
    sam2.communicate()
    
    os.system(picard + ' I=' + '{}/{}.sort.bam'.format(bam_dir, filename) +' O=' + '{}/{}.sorted.bam'.format(bam_dir, filename) + ' RGID=1 ' + ' RGPL=Illumina  RGLB=lib1 RGPU=unit1 ' + ' RGSM='+ filename )
    
    # # Index BAM
    call("{} index {}/{}.sorted.bam".format(args.samtools, bam_dir, filename).split(' '))
    
    
def consensus(args):
    """
    Using unique molecular identifiers (UMIs), duplicate reads from the same molecule are amalgamated into single-strand
    consensus sequences (SSCS). If complementary strands are present, SSCSs can be subsequently combined to form duplex
    consensus sequences (DCS).

    If 'Singleton Correction' (SC) is enabled, single reads (singletons) can be error suppressed using complementary
    strand. These corrected singletons can be merged with SSCSs to be further collapsed into DCSs + SC.

    Finally, a BAM file containing only unique molecules (i.e. no duplicates) is created by merging DCSs, remaining
    SSCSs (those that could not form DCSs), and remaining singletons (those that could not be corrected).
    """
    code_dir = os.path.dirname(os.path.realpath(__file__))

    # Change bedfile if genome is hg38
    if args.genome == 'hg38':
        # Determine code directory and set bedfile to split data
        args.bedfile = '{}/ConsensusCruncher/hg38_cytoBand.txt'.format(code_dir)

    # Create sample directory to hold consensus sequences
    identifier = os.path.basename(args.bam).split('.bam', 1)[0]
    sample_dir = '{}/{}'.format(args.c_output, identifier)

    # Check if dir exists and there's permission to write
    if not os.path.exists(sample_dir) and os.access(args.c_output, os.W_OK):
        os.makedirs(sample_dir)

    ########
    # SSCS #
    ########
    # Set variables
    os.makedirs(sample_dir + '/sscs')
    sscs = '{}/sscs/{}.sscs.bam'.format(sample_dir, identifier)
    sing = '{}/sscs/{}.singleton.bam'.format(sample_dir, identifier)

    # Run SSCS_maker
    if args.bedfile == 'False' and args.bdelim == '|':
        sscs_cmd = "{}/ConsensusCruncher/SSCS_maker.py --infile {} --outfile {} --cutoff {}".format(
            code_dir, args.bam, sscs, args.cutoff)
    elif args.bedfile == 'False' and args.bdelim != '|':
        sscs_cmd = "{}/ConsensusCruncher/SSCS_maker.py --infile {} --outfile {} --cutoff {} --bdelim {}".format(
            code_dir, args.bam, sscs, args.cutoff, args.bdelim)
    elif args.bedfile != 'False' and args.bdelim == '|':
        sscs_cmd = "{}/ConsensusCruncher/SSCS_maker.py --infile {} --outfile {} --cutoff {} --bedfile {}".format(
            code_dir, args.bam, sscs, args.cutoff, args.bedfile)
    else:
        sscs_cmd = "{}/ConsensusCruncher/SSCS_maker.py --infile {} --outfile {} --cutoff {} --bedfile {} --bdelim {}".format(
            code_dir, args.bam, sscs, args.cutoff, args.bedfile, args.bdelim)

    print(sscs_cmd)
    os.system(sscs_cmd)

    # Sort and index BAM files
    sscs = sort_index(sscs, args.samtools)
    sing = sort_index(sing, args.samtools)

    #######
    # DCS #
    #######
    # Set variables
    os.makedirs(sample_dir + '/dcs')
    dcs = '{}/dcs/{}.dcs.bam'.format(sample_dir, identifier)
    sscs_sing = '{}/dcs/{}.sscs.singleton.bam'.format(sample_dir, identifier)

    # Move stats and time tracker file to next dir
    os.rename('{}/sscs/{}.stats.txt'.format(sample_dir, identifier),
              '{}/dcs/{}.stats.txt'.format(sample_dir, identifier))
    os.rename('{}/sscs/{}.time_tracker.txt'.format(sample_dir, identifier),
              '{}/dcs/{}.time_tracker.txt'.format(sample_dir, identifier))

    # Run DCS_maker
    if args.bedfile == 'False':
        dcs_cmd = "{}/ConsensusCruncher/DCS_maker.py --infile {} --outfile {}".format(code_dir, sscs, dcs)
    else:
        dcs_cmd = "{}/ConsensusCruncher/DCS_maker.py --infile {} --outfile {} --bedfile {}".format(code_dir, sscs,
                                                                                                   dcs, args.bedfile)
    print(dcs_cmd)
    os.system(dcs_cmd)

    # Sort and index BAM files
    dcs = sort_index(dcs, args.samtools)
    sscs_sing = sort_index(sscs_sing, args.samtools)

    #############################
    # Singleton Correction (SC) #
    #############################
    if args.scorrect != 'False':
        os.makedirs(sample_dir + '/sscs_sc')
        # Move stats and time tracker file to next dir
        os.rename('{}/dcs/{}.stats.txt'.format(sample_dir, identifier),
                  '{}/sscs/{}.stats.txt'.format(sample_dir, identifier))
        os.rename('{}/dcs/{}.time_tracker.txt'.format(sample_dir, identifier),
                  '{}/sscs/{}.time_tracker.txt'.format(sample_dir, identifier))

        if args.bedfile == 'False':
            sc_cmd = "{}/ConsensusCruncher/singleton_correction.py --singleton {}".format(code_dir, sing)
        else:
            sc_cmd = "{}/ConsensusCruncher/singleton_correction.py --singleton {} --bedfile {}".format(code_dir, sing,
                                                                                                       args.bedfile)
        print(sc_cmd)
        os.system(sc_cmd)

        # Sort and index BAM files
        sscs_cor = '{}/sscs_sc/{}.sscs.correction.bam'.format(sample_dir, identifier)
        os.rename('{}/sscs/{}.sscs.correction.bam'.format(sample_dir, identifier), sscs_cor)
        sscs_cor = sort_index(sscs_cor, args.samtools)

        sing_cor = '{}/sscs_sc/{}.singleton.correction.bam'.format(sample_dir, identifier)
        os.rename('{}/sscs/{}.singleton.correction.bam'.format(sample_dir, identifier), sing_cor)
        sing_cor = sort_index(sing_cor, args.samtools)

        uncorrected = '{}/sscs_sc/{}.uncorrected.bam'.format(sample_dir, identifier)
        os.rename('{}/sscs/{}.uncorrected.bam'.format(sample_dir, identifier), uncorrected)
        uncorrected = sort_index(uncorrected, args.samtools)

        #############
        # SSCS + SC #
        #############
        # Merge corrected singletons with consensus sequences
        sscs_sc = '{}/sscs_sc/{}.sscs.sc.bam'.format(sample_dir, identifier)
        merge_sc = "{} merge {} {} {} {}".format(args.samtools, sscs_sc, sscs, sscs_cor, sing_cor)
        print(merge_sc)
        call(merge_sc.split(' '))
        sscs_sc = sort_index(sscs_sc, args.samtools)

        ############
        # DCS + SC #
        ############
        os.makedirs(sample_dir + '/dcs_sc')
        dcs_sc = '{}/dcs_sc/{}.dcs.sc.bam'.format(sample_dir, identifier)
        # Move stats and time tracker file to next dir
        os.rename('{}/sscs/{}.stats.txt'.format(sample_dir, identifier),
                  '{}/dcs_sc/{}.stats.txt'.format(sample_dir, identifier))
        os.rename('{}/sscs/{}.time_tracker.txt'.format(sample_dir, identifier),
                  '{}/dcs_sc/{}.time_tracker.txt'.format(sample_dir, identifier))

        if args.bedfile == 'False':
            dcs_sc_cmd = "{}/ConsensusCruncher/DCS_maker.py --infile {} --outfile {}".format(code_dir, sscs_sc, dcs_sc)
        else:
            dcs_sc_cmd = "{}/ConsensusCruncher/DCS_maker.py --infile {} --outfile {} --bedfile {}".format(
                code_dir, sscs_sc, dcs_sc, args.bedfile)
        print(dcs_sc_cmd)
        os.system(dcs_sc_cmd)

        # Sort and index BAM files
        dcs_sc = sort_index(dcs_sc, args.samtools)
        sscs_sc_sing = '{}/dcs_sc/{}.sscs.sc.singleton.bam'.format(sample_dir, identifier)
        sscs_sc_sing = sort_index(sscs_sc_sing, args.samtools)

        ########################
        # All Unique Molecules #
        ########################
        # Merge DCS_SC + SSCS_SC singletons + uncorrected singletons
        all_unique = '{}/dcs_sc/{}.all.unique.dcs.bam'.format(sample_dir, identifier)
        merge_all_unique = "{} merge {} {} {} {}".format(args.samtools, all_unique, dcs_sc,
                                                         sscs_sc_sing, uncorrected).split(' ')
        print(all_unique)
        call(merge_all_unique)
        all_unique = sort_index(all_unique, args.samtools)

        # Move stats and time tracker file to sample_dir
        os.rename('{}/dcs_sc/{}.stats.txt'.format(sample_dir, identifier),
                  '{}/{}.stats.txt'.format(sample_dir, identifier))
        os.rename('{}/dcs_sc/{}.time_tracker.txt'.format(sample_dir, identifier),
                  '{}/{}.time_tracker.txt'.format(sample_dir, identifier))
    else:
        os.rename('{}/dcs/{}.stats.txt'.format(sample_dir, identifier),
                  '{}/{}.stats.txt'.format(sample_dir, identifier))
        os.rename('{}/dcs/{}.time_tracker.txt'.format(sample_dir, identifier),
                  '{}/{}.time_tracker.txt'.format(sample_dir, identifier))

    # Move stats and time tracker file to sample dir
    os.rename('{}/sscs/{}_tag_fam_size.png'.format(sample_dir, identifier),
              '{}/{}_tag_fam_size.png'.format(sample_dir, identifier))
    os.rename('{}/sscs/{}.read_families.txt'.format(sample_dir, identifier),
              '{}/{}.read_families.txt'.format(sample_dir, identifier))

    # Remove intermediate files
    if args.cleanup == 'True':
        os.remove('{}/{}.time_tracker.txt'.format(sample_dir, identifier))
        os.remove('{}/sscs/{}.badReads.bam'.format(sample_dir, identifier))
        # Remove SSCSs that could not be formed into DCSs
        os.remove(sscs_sing)
        os.remove('{}/dcs/{}.sscs.singleton.sorted.bam.bai'.format(sample_dir, identifier))
        if args.scorrect != 'False':
            # Remove singleton correction files and only keep merged files
            os.remove(sing_cor)
            os.remove('{}/sscs_sc/{}.singleton.correction.sorted.bam.bai'.format(sample_dir, identifier))
            os.remove(sscs_cor)
            os.remove('{}/sscs_sc/{}.sscs.correction.sorted.bam.bai'.format(sample_dir, identifier))
            os.remove(uncorrected)
            os.remove('{}/sscs_sc/{}.uncorrected.sorted.bam.bai'.format(sample_dir, identifier))
            # Remove SSCS_SC that could not be formed into DCSs
            os.remove(sscs_sc_sing)
            os.remove('{}/dcs_sc/{}.sscs.sc.singleton.sorted.bam.bai'.format(sample_dir, identifier))


if __name__ == '__main__':
    # Set up mode parser (turn off help message, to be added later)
    main_p = argparse.ArgumentParser(add_help=False)
    main_p.add_argument('-c', '--config', default=None,
                        help="Specify config file. Commandline option overrides config file (Use config template).")

    # Parse out config file (sub_args) and other command line args (remaining_args) to override config
    sub_args, remaining_args = main_p.parse_known_args()

    # Re-initialize parser with help message enabled
    main_p = argparse.ArgumentParser(parents=[main_p], add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    sub = main_p.add_subparsers(help='sub-command help', dest='subparser_name')

    # Mode help messages
    mode_fastq2bam_help = "Extract molecular barcodes from paired-end sequencing reads using a barcode list and/or " \
                          "a barcode pattern."
    mode_consensus_help = "Almalgamate duplicate reads in BAM files into single-strand consensus sequences (SSCS) and" \
                          " duplex consensus sequences (DCS). Single reads with complementary duplex strands can also" \
                          " be corrected with 'Singleton Correction'."

    # Add subparsers
    sub_a = sub.add_parser('fastq2bam', help=mode_fastq2bam_help)
    sub_b = sub.add_parser('consensus', help=mode_consensus_help)

    # fastq2bam arg help messages
    fastq1_help = "FASTQ containing Read 1 of paired-end reads. [MANDATORY]"
    fastq2_help = "FASTQ containing Read 2 of paired-end reads. [MANDATORY]"
    output_help = "Output directory, where barcode extracted FASTQ and BAM files will be placed in " \
                  "subdirectories 'fastq_tag' and 'bamfiles' respectively (dir will be created if they " \
                  "do not exist). [MANDATORY]"
    filename_help = "Output filename. If none provided, default will extract output name by taking everything left of" \
                    " '_R'."
    bwa_help = "Path to executable bwa. [MANDATORY]"
    picard_help = "Path to executable picard add readgroups. [MANDATORY]"
    samtools_help = "Path to executable samtools. [MANDATORY]"
    ref_help = "Reference (BWA index). [MANDATORY]"
    genome_help = "Genome version (e.g. hg19 or hg38), default: hg19"
    bpattern_help = "Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases). [MANDATORY]"
    blist_help = "List of barcodes (Text file with unique barcodes on each line). [MANDATORY]"
    bdelim_help = "Delimiter before barcode in read name " \
                  "(e.g. '|' in 'HWI-D00331:196:C900FANXX:7:1110:14056:43945|TTTT')"

    # Consensus arg help messages
    bam_help = "Input BAM file with barcodes extracted into header. [MANDATORY]"
    coutput_help = "Output directory, where a folder will be created for the BAM file and consensus sequences. " \
                   "[MANDATORY]"
    scorrect_help = "Singleton correction, default: True."
    bedfile_help = "Bedfile, default: cytoBand.txt. WARNING: It is HIGHLY RECOMMENDED that you use the default " \
                   "cytoBand.txt unless you're working with genome build that is not hg19 or hg38. Then a separate " \
                   "bedfile is needed for data segmentation (file can be formatted with the bed_separator.R tool). " \
                   "For small BAM files, you may choose to turn off data splitting with '-b False' and process " \
                   "everything all at once (Division of data is only required for large data sets to offload the " \
                   "memory burden)."
    cleanup_help = "Remove intermediate files."

    # Determine code directory and set bedfile to split data
    code_dir = os.path.dirname(os.path.realpath(__file__))
    bedfile = '{}/ConsensusCruncher/hg19_cytoBand.txt'.format(code_dir)

    # Update subparsers with config
    if sub_args.config is not None:
        defaults = {"fastq1": fastq1_help,
                    "fastq2": fastq2_help,
                    "output": output_help,
                    "name": "_R",
                    "bwa": bwa_help,
                    "picard": picard_help,
                    "ref": ref_help,
                    "samtools": samtools_help,
                    "bpattern": None,
                    "blist": None,
                    "bam": bam_help,                   
                    "c_output": coutput_help,
                    "scorrect": 'True',
                    "genome": 'hg19',
                    "bedfile": bedfile,
                    "cutoff": 0.7,
                    "bdelim": '|',
                    "cleanup": cleanup_help}

        config = configparser.ConfigParser()
        config.read(sub_args.config)

        if config.has_section("fastq2bam"):
            # Add config file args to fastq2bam mode
            defaults.update(dict(config.items("fastq2bam")))
            sub_a.set_defaults(**defaults)
        if config.has_section("consensus"):
            # Add config file args to consensus mode
            defaults.update(dict(config.items("consensus")))
            sub_b.set_defaults(**defaults)

    # Parse commandline arguments
    sub_a.add_argument('--fastq1', dest='fastq1', metavar="FASTQ1", type=str, default = fastq1_help, help=fastq1_help)
    sub_a.add_argument('--fastq2', dest='fastq2', metavar="FASTQ2", type=str,default = fastq2_help,  help=fastq2_help)
    sub_a.add_argument('-o', '--output', dest='output', type=str, default = output_help,  help=output_help)
    sub_a.add_argument('-n', '--name', metavar="FILENAME", type=str, default ="_R" , help=filename_help)
    sub_a.add_argument('-b', '--bwa', metavar="BWA" , help=bwa_help, type=str)
    sub_a.add_argument('-g', '--picard', metavar="PICARD", help=picard_help, type=str)
    sub_a.add_argument('-r', '--ref', metavar="REF", help=ref_help, type=str)
    sub_a.add_argument('-s', '--samtools', metavar="SAMTOOLS",help=samtools_help, type=str)
    sub_a.add_argument('-p', '--bpattern', metavar="PATTERN",  type=str, help=bpattern_help)
    sub_a.add_argument('-l', '--blist', metavar="LIST",  type=str, help=blist_help)
    sub_a.set_defaults(func=fastq2bam)

    # Set args for 'consensus' mode
    sub_b.add_argument('-i', '--input', metavar="BAM", dest='bam', help=bam_help, type=str)
    sub_b.add_argument('-o', '--output', metavar="OUTPUT", dest='c_output', type=str, help=coutput_help)
    sub_b.add_argument('-s', '--samtools', metavar="SAMTOOLS", help=samtools_help, type=str)
    sub_b.add_argument('--scorrect', default= 'True',  help=scorrect_help, choices=['True', 'False'])
    sub_b.add_argument('-g', '--genome', metavar="VERSION", dest='genome', default = 'hg19', help=genome_help, choices=['hg19', 'hg38'])
    sub_b.add_argument('-b', '--bedfile', help=bedfile_help, type=str)
    sub_b.add_argument('--cutoff', type=float, default = 0.7, help="Consensus cut-off, default: 0.7 (70%% of reads must have the "
                                                    "same base to form a consensus).")
    sub_b.add_argument('-d', '--bdelim', metavar="DELIMITER", default ='|', type=str, help=bdelim_help)
    sub_b.add_argument('--cleanup', choices=['True', 'False'], help=cleanup_help) # Make default
    sub_b.set_defaults(func=consensus)

    # Parse args
    args = main_p.parse_args()

    # Check args
    if args.subparser_name is None:
        main_p.print_help()
    else:
        if args.subparser_name == 'fastq2bam':
            # Parse command line args to override config args
            if sub_args.config:
                sub_a.parse_known_args(remaining_args)

            # Check if required arguments provided
            if args.fastq1 is None or args.fastq2 is None or args.output is None or args.bwa is None or \
                            args.ref is None or args.samtools is None:
                sub_a.print_help()
            # Check if either barcode pattern or list is set. At least one must be provided.
            elif args.bpattern is None and args.blist is None:
                sub_a.error("At least one of -p/--bpattern or -l/--blist required.")
            # Check proper barcode design provided for barcode pattern
            elif args.bpattern is not None and re.findall(r'[^A|C|G|T|N]', args.bpattern):
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
        elif args.subparser_name == 'consensus':
            # Parse commandline args to override config args
            if sub_args.config:
                sub_b.parse_known_args(remaining_args)
            # Check if required arguments provided
            if args.bam is None or args.c_output is None or args.samtools is None:
                sub_b.print_help()
            else:
                args.func(args)
        else:
            main_p.print_help()
