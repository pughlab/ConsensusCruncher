#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

##====================================================================================
##
##  FILE:         DuplexPipeline.sh
##
##  USAGE:        DuplexPipeline.sh -i input_dir -o output_dir
##
##  OPTIONS:
##
##    -i  Input bamfile directory [MANDATORY]
##    -o  Output project directory [MANDATORY]
##    -s  Singleton correction, default: ON (use "OFF" to disable)
##    -b  Bedfile, default: cytoBand.txt
##        WARNING: It is HIGHLY RECOMMENDED that you use the default cytoBand.txt and
##        not to include your own bedfile. This option is mainly intended for non-human
##        genomes, where a separate bedfile is needed for data segmentation. If you do
##        choose to use your own bedfile, please format with the bed_separator.R tool.
##
##        For small or non-human genomes where cytobands cannot be used for segmenting the
##        data set, you may choose to turn off this option with "-b OFF" and process the
##        data all at once (Division of data is only required for large data sets to offload
##        the memory burden).
##
##    -c  Consensus cut-off, default: 0.7 (70% of reads must have the same base to form
##        a consensus)
##    -q  qusb directory, default: output/qsub
##    -h  Show this message
##
##  DESCRIPTION:
##
##  This script amalgamates duplicate reads in bamfiles into single-strand consensus
##  sequences (SSCS), which are subsequently combined into duplex consensus sequences
##  (DCS). Singletons (reads lacking duplicate sequences) are corrected, combined
##  with SSCS to form SSCS + SC, and further collapsed to form DCS + SC. Finally,
##  files containing all unique molecules (a.k.a. no duplicates) are created for SSCS
##  and DCS.
##
##  Note: Script will create a "consensus" directory under the project directory and
##  sub-directories corresponding to each bamfile in the input directory.
##
##  WARNING: Please change qsub parameters according to your cluster commands; default
##  qsub -q highmem.q SCRIPT
##
##====================================================================================

usage()
{
cat << EOF

  FILE:         DuplexPipeline.sh

  USAGE:        DuplexPipeline.sh -i input_dir -o output_dir

  OPTIONS:

    -i  Input bamfile directory [MANDATORY]
    -o  Output project directory [MANDATORY]
    -s  Singleton correction, default: ON (use "OFF" to disable)
    -b  Bedfile, default: cytoBand.txt
        WARNING: It is HIGHLY RECOMMENDED that you use the default cytoBand.txt and
        not to include your own bedfile. This option is mainly intended for non-human
        genomes, where a separate bedfile is needed for data segmentation. If you do
        choose to use your own bedfile, please format with the bed_separator.R tool.

        For small or non-human genomes where cytobands cannot be used for segmenting the
        data set, you may choose to turn off this option with "-b OFF" and process the
        data all at once (Division of data is only required for large data sets to offload
        the memory burden).

    -c  Consensus cut-off, default: 0.7 (70% of reads must have the same base to form
        a consensus)
    -q  qusb directory, default: output/qsub
    -h  Show this message

  DESCRIPTION:

  This script amalgamates duplicate reads in bamfiles into single-strand consensus
  sequences (SSCS), which are subsequently combined into duplex consensus sequences
  (DCS). Singletons (reads lacking duplicate sequences) are corrected, combined
  with SSCS to form SSCS + SC, and further collapsed to form DCS + SC. Finally,
  files containing all unique molecules (a.k.a. no duplicates) are created for SSCS
  and DCS.

  Note: Script will create a "consensus" directory under the project directory and
  sub-directories corresponding to each bamfile in the input directory.

  WARNING: Please change qsub parameters according to your cluster commands; default
  qsub -q highmem.q SCRIPT

EOF
}

################
#    Set-up    #
################
while getopts "hi:o:s:b:c:q:" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             INPUT=$OPTARG
             ;;
         o)
             OUTPUT=$OPTARG
             ;;
         s)
             SINGCOR=$OPTARG
             ;;
         b)
             BEDFILE=$OPTARG
             ;;
         c)
             CUTOFF=$OPTARG
             ;;
         q)
             QSUBDIR=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $INPUT ]] || [[ -z $OUTPUT ]]
then
     usage
     exit 1
fi


################
#    Modules   #
################
#load modules
module load python3/3.4.3
module load pysam
module load samtools/1.2
module load picard/2.4.1
module load java/8 

##check necessary modules have been loaded
tools="python pysam samtools picard java"
modules=$(module list 2>&1)
for i in $tools
do
    loaded="$(echo $modules|grep -oP "$i.*?([[:space:]]|$)")"
    n=$(wc -l <<<"$loaded")
    if [ "$n" -gt "1" ]; then
        >&2 echo "Error: multiple $i modules were loaded. I am confused."
        exit 1
    fi

    if [[ -z $loaded ]]
    then
        echo Please load module for $i
        exit 1
    else
        eval "${i}_version=$loaded"
        echo "I am going to use $loaded"
    fi
done

# Set script directory
code_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
code_dir=$code_dir/consensus_scripts

###############
#  Variables  #
###############
# Set bedfile variable
if [[ -z $BEDFILE ]]; then
    BEDFILE=$code_dir/cytoBand.txt
elif [[ $BEDFILE == "OFF" ]]; then
    unset BEDFILE
fi

# Set Cut-off if none provided
if [[ -z $CUTOFF ]]; then
    CUTOFF=0.7
fi


################
#  Create dir  #
################
: ${OUTPUT:=$OUTPUT/consensus}
[ -d $OUTPUT ] || mkdir -p $OUTPUT
: ${QSUBDIR:=$OUTPUT/qsub}
[ -d $QSUBDIR ] || mkdir $QSUBDIR


####################
#  Duplex Pipeline #
####################
for bamfile in $(ls $INPUT | grep -v bai); do
    ###############
    #  Sample ID  #
    ###############
    identifier=${bamfile//\.bam/}

    # Create sample directory
    : ${SAMPDIR:=$OUTPUT/$identifier}
    [ -d $SAMPDIR ] || mkdir $SAMPDIR

    # Set-up sh script
    echo "#!/bin/bash" > $QSUBDIR/$identifier.sh
    echo -e "#$ -S /bin/bash\n#$ -cwd\n#$ -N\nmodule load $python_version\nmodule load $pysam_version\nmodule load $samtools_version\nmodule load $picard_version\nmodule load $java_version\n">>$QSUBDIR/$identifier.sh

    echo -e "cd $SAMPDIR\n" >> $QSUBDIR/$identifier.sh

    ################
    #     SSCS     #
    ################
    if [[ -z $BEDFILE ]]; then
        echo -e "python3 $code_dir/SSCS_maker.py  --cutoff $CUTOFF --infile $INPUT/$bamfile --outfile $SAMPDIR/$identifier.sscs.bam" >> $QSUBDIR/$identifier.sh
    else
        echo -e "python3 $code_dir/SSCS_maker.py  --cutoff $CUTOFF --infile $INPUT/$bamfile --outfile $SAMPDIR/$identifier.sscs.bam --bedfile $BEDFILE" >> $QSUBDIR/$identifier.sh
    fi
    # sort and index SSCS
    echo -e "samtools view -bu $identifier.sscs.bam | samtools sort - $identifier.sscs.sorted" >> $QSUBDIR/$identifier.sh
    echo -e "samtools index $identifier.sscs.sorted.bam" >> $QSUBDIR/$identifier.sh
    echo -e "rm $identifier.sscs.bam" >> $QSUBDIR/$identifier.sh
    # sort and index Singletons
    echo -e "samtools view -bu $identifier.singleton.bam | samtools sort - $identifier.singleton.sorted" >> $QSUBDIR/$identifier.sh
    echo -e "samtools index $identifier.singleton.sorted.bam" >> $QSUBDIR/$identifier.sh
    echo -e "rm $identifier.singleton.bam" >> $QSUBDIR/$identifier.sh

    ###############
    #     DCS     #
    ###############
    if [[ -z $BEDFILE ]]; then
        echo -e "python3 $code_dir/DCS_maker.py --infile $identifier.sscs.sorted.bam --outfile $identifier.dcs.bam" >> $QSUBDIR/$identifier.sh
    else
        echo -e "python3 $code_dir/DCS_maker.py --infile $identifier.sscs.sorted.bam --outfile $identifier.dcs.bam --bedfile $BEDFILE" >> $QSUBDIR/$identifier.sh
    fi
    # sort and index DCS
    echo -e "samtools view -bu $identifier.dcs.bam | samtools sort - $identifier.dcs.sorted" >> $QSUBDIR/$identifier.sh
    echo -e "samtools index $identifier.dcs.sorted.bam" >> $QSUBDIR/$identifier.sh
    echo -e "rm $identifier.dcs.bam" >> $QSUBDIR/$identifier.sh
    # sort and index SSCS singletons
    echo -e "samtools view -bu $identifier.sscs.singleton.bam | samtools sort - $identifier.sscs.singleton.sorted" >> $QSUBDIR/$identifier.sh
    echo -e "samtools index $identifier.sscs.singleton.sorted.bam" >> $QSUBDIR/$identifier.sh
    echo -e "rm $identifier.sscs.singleton.bam" >> $QSUBDIR/$identifier.sh

    #############################
    # Singleton Correction (SC) #
    #############################
    if [[ -z $SINGCOR ]] || [[ $SINGCOR == "ON" ]]; then
        if [[ -z $BEDFILE ]]; then
            echo -e "python3 $code_dir/singleton_correction.py --singleton $identifier.singleton.sorted.bam" >> $QSUBDIR/$identifier.sh
        else
            echo -e "python3 $code_dir/singleton_correction.py --singleton $identifier.singleton.sorted.bam --bedfile $BEDFILE" >> $QSUBDIR/$identifier.sh
        fi
        # sort and index Singletons Correction by SSCS
        echo -e "samtools view -bu $identifier.sscs.rescue.bam | samtools sort - $identifier.sscs.rescue.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.sscs.rescue.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.sscs.rescue.bam" >> $QSUBDIR/$identifier.sh
        # sort and index Singletons Correction by Singletons
        echo -e "samtools view -bu $identifier.singleton.rescue.bam | samtools sort - $identifier.singleton.rescue.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.singleton.rescue.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.singleton.rescue.bam" >> $QSUBDIR/$identifier.sh
        # sort and index remaining singletons
        echo -e "samtools view -bu $identifier.rescue.remaining.bam | samtools sort - $identifier.rescue.remaining.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.rescue.remaining.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.rescue.remaining.bam" >> $QSUBDIR/$identifier.sh

        #############
        # SSCS + SC #
        #############
        echo -e "java -jar $picard_dir/picard.jar MergeSamFiles I=$identifier.sscs.sorted.bam I=$identifier.sscs.rescue.sorted.bam I=$identifier.singleton.rescue.sorted.bam O=$identifier.sscs.sc.bam" >> $QSUBDIR/$identifier.sh
        echo -e "samtools view -bu $identifier.sscs.sc.bam | samtools sort - $identifier.sscs.sc.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.sscs.sc.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.sscs.sc.bam" >> $QSUBDIR/$identifier.sh

        ############
        # DCS + SC #
        ############
        # Construct DCS from SSCS + SC
        if [[ -z $BEDFILE ]]; then
            echo -e "python3 $code_dir/DCS_maker.py --infile $identifier.sscs.sc.sorted.bam --outfile $identifier.dcs.sc.bam" >> $QSUBDIR/$identifier.sh
        else
            echo -e "python3 $code_dir/DCS_maker.py --infile $identifier.sscs.sc.sorted.bam --outfile $identifier.dcs.sc.bam --bedfile $BEDFILE" >> $QSUBDIR/$identifier.sh
        fi
        # sort and index DCS SC
        echo -e "samtools view -bu $identifier.dcs.sc.bam | samtools sort - $identifier.dcs.sc.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.dcs.sc.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.dcs.sc.bam" >> $QSUBDIR/$identifier.sh
        # sort and index SSCS SC Singletons (a.k.a. remaining SSCS + SC that could not be made into DCS)
        echo -e "samtools view -bu $identifier.sscs.sc.singleton.bam | samtools sort - $identifier.sscs.sc.singleton.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.sscs.sc.singleton.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.sscs.sc.singleton.bam" >> $QSUBDIR/$identifier.sh

        ########################
        # All Unique Molecules #
        ########################
        # === SSCS + SC + remaining singletons ===
        echo -e "java -jar $picard_dir/picard.jar MergeSamFiles I=$identifier.sscs.sorted.bam I=$identifier.sscs.rescue.sorted.bam I=$identifier.singleton.rescue.sorted.bam I=$identifier.rescue.remaining.sorted.bam O=$identifier.all.unique.sscs.bam" >> $QSUBDIR/$identifier.sh
        # sort and index all unique SSCS molecules
        echo -e "samtools view -bu $identifier.all.unique.sscs.bam | samtools sort - $identifier.all.unique.sscs.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.all.unique.sscs.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.all.unique.sscs.bam" >> $QSUBDIR/$identifier.sh

        # === DCS (SSCS_SC) + SSCS SC singletons + remaining singletons ===
        echo -e "java -jar $picard_dir/picard.jar MergeSamFiles I=$identifier.dcs.sc.sorted.bam I=$identifier.sscs.sc.singleton.sorted.bam I=$identifier.rescue.remaining.sorted.bam O=$identifier.all.unique.dcs.bam" >> $QSUBDIR/$identifier.sh
        # sort and index all unique DCS molecules
        echo -e "samtools view -bu $identifier.all.unique.dcs.bam | samtools sort - $identifier.all.unique.dcs.sorted" >> $QSUBDIR/$identifier.sh
        echo -e "samtools index $identifier.all.unique.dcs.sorted.bam" >> $QSUBDIR/$identifier.sh
        echo -e "rm $identifier.all.unique.dcs.bam" >> $QSUBDIR/$identifier.sh

        #####################
        # Organize SC files #
        #####################
        echo -e "mkdir sscs_SC" >> $QSUBDIR/$identifier.sh
        echo -e "mkdir dcs_SC" >> $QSUBDIR/$identifier.sh

        # DCS + SC
        echo -e "mv $identifier.all.unique.dcs.sorted* dcs_SC" >> $QSUBDIR/$identifier.sh
        echo -e "mv $identifier.dcs.sc* dcs_SC" >> $QSUBDIR/$identifier.sh
        echo -e "mv $identifier.sscs.sc.singleton* dcs_SC" >> $QSUBDIR/$identifier.sh

        # SSCS + SC
        echo -e "mv $identifier.all.unique.sscs* sscs_SC" >> $QSUBDIR/$identifier.sh
        echo -e "mv $identifier.sscs.sc* sscs_SC" >> $QSUBDIR/$identifier.sh
        echo -e "mv $identifier.*rescue* sscs_SC" >> $QSUBDIR/$identifier.sh

    fi

    #########################
    # Organize files by dir #
    #########################
    echo -e "mkdir sscs" >> $QSUBDIR/$identifier.sh
    echo -e "mkdir dcs" >> $QSUBDIR/$identifier.sh


    # DCS
    echo -e "mv $identifier.dcs.sorted.* dcs" >> $QSUBDIR/$identifier.sh
    echo -e "mv $identifier.sscs.singleton.sorted.* dcs" >> $QSUBDIR/$identifier.sh

    # SSCS
    echo -e "mv $identifier* sscs" >> $QSUBDIR/$identifier.sh

    # Keep stats and family size plot in sample directory
    echo -e "mv ./sscs/*stats.txt ." >> $QSUBDIR/$identifier.sh
    echo -e "mv ./sscs/*png ." >> $QSUBDIR/$identifier.sh

    ########
    # QSUB #
    ########
    cd $QSUBDIR
    qsub -q highmem.q $QSUBDIR/$identifier.sh

done