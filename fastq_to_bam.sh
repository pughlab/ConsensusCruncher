#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

##====================================================================================
##
##  FILE:         fastq_to_bam.sh
##
##  USAGE:        fastq_to_bam.sh -i fastq -o bwa_out -p project -t 2 -s 1 -f T
##
##  OPTIONS:
##
##    -i  Input directory.[MANDATORY]
##    -o  Output directory [MANDATORY]
##    -p  Project name [MANDATORY]
##    -t  Tag length [MANDATORY]
##    -s  Spacer length [MANDATORY]
##    -f  Spacer Filter (e.g. "T" will filter out spacers that are non-T)
##    -q  qusb directory, default: output/qsub
##    -h  Show this message
##
##  DESCRIPTION:
##
##  This script removes molecular barcode tags  and spacers from unzipped FASTQ files,
##  which are found by searching the source directory for files with "*R1". Tag removed
##  FASTQ files are subsequently aligned with BWA mem.
##
##====================================================================================

usage()
{
cat << EOF

  FILE:         fastq_to_bam.sh

  USAGE:        fastq_to_bam.sh -i fastq -o bwa_out -p project -t 2 -s 1 -f T

  OPTIONS:

    -i  Source directory.[MANDATORY]
    -o  Output directory [MANDATORY]
    -p  Project name [MANDATORY]
    -t  Tag length [MANDATORY]
    -s  Spacer length [MANDATORY]
    -f  Spacer filter (e.g. "T" will filter out spacers that are non-T)
    -q  qusb directory, default: output/qsub
    -h  Show this message

  DESCRIPTION:

  This script removes molecular barcode tags and spacers from unzipped FASTQ files,
  which are found by searching the source directory for files with "*R1". Tag removed
  FASTQ files are subsequently aligned with BWA mem.

EOF
}

################
#    Set-up    #
################

while getopts "hi:o:p:t:s:f:q:" OPTION
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
         p)
             PROJECT=$OPTARG
             ;;
         t)
             TAGLEN=$OPTARG
             ;;
         s)
             SPACERLEN=$OPTARG
             ;;
         f)
             SPACERFILT=$OPTARG
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

if [[ -z $INPUT ]] || [[ -z $OUTPUT ]] || [[ -z $PROJECT ]] || [[ -z $TAGLEN ]] || [[ -z $SPACERLEN ]]
then
     usage
     exit 1
fi


################
#    Modules   #
################
## Load modules
module load python
module load samtools
module load bwa
BWAINDEX="/mnt/work1/data/genomes/human/hg19/iGenomes/Sequence/BWAIndex/genome.fa"

##check necessary modules have been loaded
tools="python bwa samtools"
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

##check necessary variables have been set
variables="BWAINDEX"
for i in $variables; do
    if [ -z ${!i} ]; then
        echo "Please set variable $i"
        exit 1
    else
        echo "I am going to use $i=${!i}"
    fi
done

# Set code directory
code_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


################
#  Create dir  #
################
[ -d $OUTPUT ] || mkdir -p $OUTPUT
: ${QSUBDIR:=$OUTPUT/qsub}
[ -d $QSUBDIR ] || mkdir $QSUBDIR
: ${TAGDIR:=$OUTPUT/fastq_tag}
[ -d $TAGDIR ] || mkdir $TAGDIR


#####################
# Create sh scripts #
#####################

for R1_file in $( ls $INPUT | grep R1); do
    ################
    #  Set-up IDs  #
    ################
    R2_file=${R1_file//R1/R2}
    filename=${R1_file//_R1/}
    filename=${filename//.fastq/}

    IFS='_' read -a fname_array<<<"$filename"  # Separate filename by "_"
    lane="$(printf -- "%s\n" "${fname_array[@]}" | grep L00)"
    lane_i="$(echo ${fname_array[@]/$lane/-} | cut -d- -f1 | wc -w | tr -d ' ')"
    barcode_i="$(expr $lane_i - 1)"
    barcode=${fname_array[$barcode_i]}

    echo -e "#/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n\nmodule load $python_version\nmodule load $bwa_version\nmodule load $samtools_version\n" > $QSUBDIR/$filename.sh

    #################
    #  Unzip files  #
    #################
    if [[ $R1_file = *"gz" ]]; then
        # Make unzipped file dir
        : ${UNZIPDIR:=$OUTPUT/fastq_unzip}
        [ -d $UNZIPDIR ] || mkdir $UNZIPDIR

        # Create tmp unzipped file
        R1_unzip=${R1_file//.gz/}
        echo -e "zcat $INPUT/$R1_file > $UNZIPDIR/$R1_unzip" >> $QSUBDIR/$filename.sh
        R2_unzip=${R2_file//.gz/}
        echo -e "zcat $INPUT/$R2_file > $UNZIPDIR/$R2_unzip" >> $QSUBDIR/$filename.sh

        R1="$UNZIPDIR/$R1_unzip"
        R2="$UNZIPDIR/$R2_unzip"
    else
        R1="$INPUT/$R1_file"
        R2="$INPUT/$R2_file"
    fi

    #################
    #  Remove tags  #
    #################
    # Change directory so tag stats file will be created in $TAGDIR
    cd $TAGDIR

    # Check if there's a spacer filter
    if [ -z $SPACERFILT ]; then
        echo -e "python $code_dir/consensus_scripts/tag_to_header.py --infile1 $R1 --infile2 $R2 --outprefix $TAGDIR/$filename --taglen $TAGLEN --spacerlen $SPACERLEN --filtspacer $SPACERFILT \n" >> $QSUBDIR/$filename.sh
    else
        echo -e "python $code_dir/consensus_scripts/tag_to_header.py --infile1 $R1 --infile2 $R2 --outprefix $TAGDIR/$filename --taglen $TAGLEN --spacerlen $SPACERLEN \n" >> $QSUBDIR/$filename.sh
    fi

    #################
    #  Align reads  #
    #################
    echo -e "bwa mem -M -t4 -R '@RG\tID:1\tSM:$filename\tPL:Illumina\tPU:$barcode.$lane\tLB:$PROJECT' $BWAINDEX $TAGDIR/$filename.seq1.smi.fq $TAGDIR/$filename.seq2.smi.fq > $OUTPUT/$filename.sam \n" >>$QSUBDIR/$filename.sh

    # Convert to BAM format and sort by positions
    echo -e "samtools view -bhS $OUTPUT/$filename.sam | samtools sort -@4 - $OUTPUT/$filename \n" >> $QSUBDIR/$filename.sh
    echo -e "samtools index $OUTPUT/$filename.bam" >> $QSUBDIR/$filename.sh

    # Remove sam file
    echo -e "rm $OUTPUT/$filename.sam" >> $QSUBDIR/$filename.sh

    #########################
    # Remove unzipped files #
    #########################
    # Don't want to remove entire folder as there may be other files used by other scripts
    if [ -z $UNZIPDIR ]; then
        echo -e "rm $UNZIPDIR/$R1_unzip" >> $QSUBDIR/$filename.sh
        echo -e "rm $UNZIPDIR/$R2_unzip" >> $QSUBDIR/$filename.sh
    fi

    cd $QSUBDIR
    qsub $QSUBDIR/$filename.sh
done
