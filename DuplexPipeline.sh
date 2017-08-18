#/bin/bash
#$ -S /bin/bash
#$ -cwd

################
#    SET-UP    #
################
sample_dir=$1
input_bam=$2
bedfile=$3
codedir=$4

################
# Load Modules #
################
module load samtools/1.2
module load python3/3.4.3
module load picard/2.4.1 
module load java/8 

###############
#  Sample ID  #
###############
identifier=${input_bam//\.bam/}
identifier=${identifier##*/}
echo $identifier

cd $sample_dir

################
#     SSCS     #
################
python3 $codedir/SSCS_maker.py  --cutoff 0.7 --infile $input_bam --outfile $identifier.sscs.bam --bedfile $bedfile
# sort and index SSCS 
samtools view -bu $identifier.sscs.bam | samtools sort - $identifier.sscs.sorted
samtools index $identifier.sscs.sorted.bam 
rm $identifier.sscs.bam 
# sort and index Singletons
samtools view -bu $identifier.singleton.bam | samtools sort - $identifier.singleton.sorted 
samtools index $identifier.singleton.sorted.bam
rm $identifier.singleton.bam

###############
#     DCS     #
###############
python3 $codedir/DCS_maker.py --infile $identifier.sscs.sorted.bam --outfile $identifier.dcs.bam --bedfile $bedfile
# sort and index DCS 
samtools view -bu $identifier.dcs.bam | samtools sort - $identifier.dcs.sorted
samtools index $identifier.dcs.sorted.bam
rm $identifier.dcs.bam
# sort and index SSCS singletons
samtools view -bu $identifier.sscs.singleton.bam | samtools sort - $identifier.sscs.singleton.sorted
samtools index $identifier.sscs.singleton.sorted.bam
rm $identifier.sscs.singleton.bam

#########################
# Singleton Rescue (SR) #
#########################
python3 $codedir/singleton_strand_rescue.py --singleton $identifier.singleton.sorted.bam --bedfile $bedfile
# sort and index Singletons rescued by SSCS
samtools view -bu $identifier.sscs.rescue.bam | samtools sort - $identifier.sscs.rescue.sorted
samtools index $identifier.sscs.rescue.sorted.bam
rm $identifier.sscs.rescue.bam 
# sort and index Singletons rescued by Singletons
samtools view -bu $identifier.singleton.rescue.bam | samtools sort - $identifier.singleton.rescue.sorted
samtools index $identifier.singleton.rescue.sorted.bam
rm $identifier.singleton.rescue.bam
# sort and index remaining singletons 
samtools view -bu $identifier.rescue.remaining.bam | samtools sort - $identifier.rescue.remaining.sorted
samtools index $identifier.rescue.remaining.sorted.bam
rm $identifier.rescue.remaining.bam

#############
# SSCS + SR #
#############
java -jar $picard_dir/picard.jar MergeSamFiles I=$identifier.sscs.sorted.bam \
											   I=$identifier.sscs.rescue.sorted.bam \
  											   I=$identifier.singleton.rescue.sorted.bam \
											   O=$identifier.sscs.sr.bam

samtools view -bu $identifier.sscs.sr.bam | samtools sort - $identifier.sscs.sr.sorted
samtools index $identifier.sscs.sr.sorted.bam
rm $identifier.sscs.sr.bam

######################
# DCS from SSCS + SR #
######################
python3 $codedir/DCS_maker.py --infile $identifier.sscs.sr.sorted.bam --outfile $identifier.dcs.sr.bam --bedfile $bedfile
# sort and index dcs sr
samtools view -bu $identifier.dcs.sr.bam | samtools sort - $identifier.dcs.sr.sorted
samtools index $identifier.dcs.sr.sorted.bam 
rm $identifier.dcs.sr.bam
# sort and index sscs sr singletons
samtools view -bu $identifier.sscs.sr.singleton.bam | samtools sort - $identifier.sscs.sr.singleton.sorted
samtools index $identifier.sscs.sr.singleton.sorted.bam 
rm $identifier.sscs.sr.singleton.bam

########################
# All Unique Molecules #
########################
# === SSCS + SR + remaining singletons ===
java -jar $picard_dir/picard.jar MergeSamFiles I=$identifier.sscs.sorted.bam \
											   I=$identifier.sscs.rescue.sorted.bam \
											   I=$identifier.singleton.rescue.sorted.bam \
											   I=$identifier.rescue.remaining.sorted.bam \
											   O=$identifier.all.unique.sscs.bam
# sort and index all unique SSCS molecules
samtools view -bu $identifier.all.unique.sscs.bam | samtools sort - $identifier.all.unique.sscs.sorted
samtools index $identifier.all.unique.sscs.sorted.bam
rm $identifier.all.unique.sscs.bam
		
		
# === DCS (SSCS_SR) + SSCS SR singletons + remaining singletons ===
java -jar $picard_dir/picard.jar MergeSamFiles I=$identifier.dcs.sr.sorted.bam \
											   I=$identifier.sscs.sr.singleton.sorted.bam \
											   I=$identifier.rescue.remaining.sorted.bam \
											   O=$identifier.all.unique.dcs.bam
# sort and index all unique DCS molecules
samtools view -bu $identifier.all.unique.dcs.bam | samtools sort - $identifier.all.unique.dcs.sorted
samtools index $identifier.all.unique.dcs.sorted.bam
rm $identifier.all.unique.dcs.bam

#########################
# Organize files by dir #
#########################
mkdir sscs
mkdir sscs_SR
mkdir dcs
mkdir dcs_SR

# DCS SR
mv $identifier.all.unique.dcs.sorted* dcs_SR
mv $identifier.dcs.sr* dcs_SR
mv $identifier.sscs.sr.singleton* dcs_SR

# DCS
mv $identifier.dcs.sorted.* dcs
mv $identifier.sscs.singleton.sorted.* dcs

# SSCS SR
mv $identifier.all.unique.sscs* sscs_SR
mv $identifier.sscs.sr* sscs_SR
mv $identifier.*rescue* sscs_SR

# SSCS
mv $identifier* sscs

# Keep stats and family size plot in main directory
mv ./sscs/*stats.txt .
mv ./sscs/*png .

