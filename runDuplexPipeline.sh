#/bin/bash
#$ -S /bin/bash
#$ -cwd

################
#    SET-UP    #
################
project_dir=$1
cwd=$project_dir
bamdir=$2
codedir='/mnt/work1/users/pughlab/projects/Duplext_sequencing/nwang_scripts/pl_duplex_sequencing'
if [ $# -eq 3 ]; then
    bedfile=$3
else
	bedfile=$codedir/cytoBand.txt
fi

cd $project_dir
mkdir consensus

#######################
# run Duplex Pipeline #
#######################
cd $bamdir
for bamfile in $(ls *.bam); do
	qsub -q highmem.q $codedir/DuplexScriptGen.sh $project_dir $bamdir/$bamfile $bedfile $codedir/consensus_scripts
done
