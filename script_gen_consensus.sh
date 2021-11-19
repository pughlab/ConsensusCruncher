INPUT=#IndelRealigned bams
OUTPUT=# Output directory
CODEDIR="/mnt/work1/users/pughlab/bin/ConsensusCruncher/" # ConsensusCruncher code directory
CONFIG="/mnt/work1/users/pughlab/projects/WALDENSTORM/ConsensusCruncher/config.ini"  # Path to config file
CYTOBAND=$CODEDIR/ConsensusCruncher/hg38_cytoBand.txt # Textfile to separate

QSUBDIR=$OUTPUT/consensus/qsub

mkdir $OUTPUT/consensus
mkdir $QSUBDIR
for filename in $(ls $INPUT | grep -v bai); do
echo -e "#!/bin/bash\n \nmodule load python3/3.7.5\n" > $QSUBDIR/$filename.sh
###############
    #  consensus  #
    ###############
    echo -e "python3 $CODEDIR/ConsensusCruncher.py -c $CONFIG consensus -i $INPUT/$filename -o $OUTPUT/consensus -b $CYTOBAND" >> $QSUBDIR/$filename.sh
