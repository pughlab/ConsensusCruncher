INPUT=#Input directory where fastqs are sitting
OUTPUT=# output directory
BLIST=/cluster/projects/pughlab/references/UMI_list/IDT_dual_index.txt  # comment out this line if you are using Bratman barcodes
CODEDIR=##ConsensusCruncher Code directory 
PATTERN="NNT" ## comment out this line if you are using IDT barcoces
CONFIG=config.ini
CYTOBAND=$CODEDIR/ConsensusCruncher/hg38_cytoBand.txt  # Textfile to separate
QSUBDIR=$OUTPUT/consensus/qsub

mkdir $OUTPUT/consensus
mkdir $QSUBDIR

for R1_file in $( ls $INPUT | grep R1); do
 
     ################
    #  Set-up IDs  #
    ################
    R2_file=${R1_file//R1/R2}
    filename="$(echo $R1_file | sed 's/_R.*//')"  # Note file naming should be tailored as it currently removes everything right of _"R" e.g. sample_R1.fastq -> sample.fastq
    echo -e "#!/bin/bash\n#$ -S /bin/bash\n#$ -cwd\n\nmodule load python3/3.7.2\n" > $QSUBDIR/$filename.sh
    echo -e "\nmodule load picard/2.10.9\nmodule load igenome-human/hg38\n" >> $QSUBDIR/$filename.sh

    ###############
    #  fastq2bam  #
    ###############
    echo -e "python3 $CODEDIR/ConsensusCruncher.py  -c $CONFIG fastq2bam --fastq1 $INPUT/$R1_file --fastq2 $INPUT/$R2_file --output $OUTPUT -p $PATTERN" >> $QSUBDIR/$filename.sh
done
