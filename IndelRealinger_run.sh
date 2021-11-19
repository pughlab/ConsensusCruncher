#!/bin/sh
bam=#bamfile name ( you can loop through a list of bam file by doing for bam in $(cat bamfiles);do 
module load gatk/3.3-0 
module load igenome-human/hg38
input=## bam directory from the fastqTobam steph
bed=## bed file
cwd=## currecnt working directory where your processed bam files will go to
sites=/cluster/projects/pughlab/references/broad_homo_sapien_known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz

mkdir -p $cwd/intervals
mkdir -p $cwd/output
gatk -T  RealignerTargetCreator \
    -R $REF  \
    -L $bed \
    -known $sites \
    -I  $input/$1  \
    -o $cwd/intervals/${1}.bed
   
java -jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner \
    -R $REF \
    -targetIntervals $cwd/intervals/${1}.bed  \
    -known $sites -I /cluster/projects/pughlab/projects/M4/ConsensusCruncher/processed/all_bams/$bam \
    -o $output/$1
