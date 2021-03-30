#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

# Generate ConsensusCruncher sh scripts for each bam file to qsub (Please load appropriate python module)

INPUT=$1  # Bamfile directory
OUTPUT=$2  # Output directory
PATTERN=$3  # Barcode pattern, e.g. NNT or NNGCT
CODEDIR= # ConsensusCruncher code directory
CONFIG=  # Path to config file
CYTOBAND=$CODEDIR/ConsensusCruncher/hg19_cytoBand.txt # Textfile to separate

QSUBDIR=$OUTPUT/consensus/qsub

mkdir $OUTPUT/consensus
mkdir $QSUBDIR

for R1_file in $( ls $INPUT | grep R1); do
    ################
    #  Set-up IDs  #
    ################
    R2_file=${R1_file//R1/R2}
    filename="$(echo $R1_file | sed 's/_R.*//')"  # Note file naming should be tailored as it currently removes everything right of _"R" e.g. sample_R1.fastq -> sample.fastq

    echo -e "#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n\nmodule load python3/3.4.3\n" > $QSUBDIR/$filename.sh
    echo -e "\nmodule load picard/2.10.9 \n" >> $QSUBDIR/$filename.sh
    ###############
    #  fastq2bam  #
    ###############
    echo -e "python3 $CODEDIR/ConsensusCruncher.py -c $CONFIG fastq2bam --fastq1 $INPUT/$R1_file --fastq2 $INPUT/$R2_file --output $OUTPUT -p $PATTERN" >> $QSUBDIR/$filename.sh

    ###############
    #  consensus  #
    ###############
    echo -e "python3 $CODEDIR/ConsensusCruncher.py -c $CONFIG consensus -i $OUTPUT/bamfiles/$filename.sorted.bam -o $OUTPUT/consensus -b $CYTOBAND" >> $QSUBDIR/$filename.sh

done
