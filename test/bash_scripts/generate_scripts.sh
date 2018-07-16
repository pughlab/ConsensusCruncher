#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

# Generate ConsensusCruncher sh scripts for each bam file to qsub (Please load appropriate python module)

INPUT=$1  # Directory containing FASTQ files
OUTPUT=$2  # Output project directory, new folders and files will be created here
QSUBDIR=$OUTPUT/qsub
code_dir='PLEASE INSERT PATH TO GIT DIRECTORY' 

mkdir $QSUBDIR
mkdir $OUTPUT'/consensus'

for R1_file in $( ls $INPUT | grep R1); do
    ################
    #  Set-up IDs  #
    ################
    R2_file=${R1_file//R1/R2}
    filename="$(echo $R1_file | sed 's/_R.*//')"  # Note file naming should be tailored as it currently removes everything right of _"R" e.g. sample_R1.fastq -> sample.fastq

    echo -e "#/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n\nmodule load python3/3.4.3\n" > $QSUBDIR/$filename.sh

    ###############
    #  fastq2bam  #
    ###############
    echo -e "python3 $code_dir/ConsensusCruncher.py -c $code_dir/config.ini fastq2bam --fastq1 $INPUT/$R1_file --fastq2 $INPUT/$R2_file --output $OUTPUT" >> $QSUBDIR/$filename.sh

    ###############
    #  consensus  #
    ###############
    echo -e "python3 $code_dir/ConsensusCruncher.py -c $code_dir/config.ini consensus -i $OUTPUT/bamfiles/$filename.bam" >> $QSUBDIR/$filename.sh

done

