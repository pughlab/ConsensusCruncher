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
	# WARNING: Please format bedfile using bed_separator.R tool before using for duplex pipeline. 
	# It is HIGHLY RECOMMENDED to use the default cytoBand.txt and not to include your own bedfile.
    bedfile=$3  
else
	bedfile=$codedir/consensus_scripts/cytoBand.txt
fi

cd $project_dir
mkdir consensus

#######################
# run Duplex Pipeline #
#######################
cd $bamdir
for bamfile in $(ls *bam); do
	identifier=${bamfile//\.bam/}
	echo $identifier

	cd $project_dir/consensus
	mkdir 'Duplex_'$identifier
	cd 'Duplex_'$identifier
	sample_dir=$project_dir/consensus/'Duplex_'$identifier

	qsub -q highmem.q $codedir/DuplexPipeline.sh $sample_dir $bamdir/$bamfile $bedfile $codedir/consensus_scripts
done
