#!/bin/sh
module load gatk/3.3-0 
module load igenome-human/hg38
input=/cluster/projects/pughlab/projects/M4/ConsensusCruncher/processed/Genes_bams
bed=/cluster/projects/pughlab/projects/M4/Bed_files/pughlab_panels/MYLSTONE2_Gene_GRCh38.bed
cwd=/cluster/projects/pughlab/projects/M4/ConsensusCruncher/processed
sites=/cluster/projects/pughlab/references/broad_homo_sapien_known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz

mkdir -p $cwd/intervals

gatk -T  RealignerTargetCreator \
    -R $REF  \
    -L $bed \
    -known $sites \
    -I  $input/$1  \
    -o $cwd/intervals/${1}.bed
   
java -jar $gatk_dir/GenomeAnalysisTK.jar -T IndelRealigner \
    -R $REF \
    -targetIntervals $cwd/intervals/${1}.bed  \
    -known $sites -I /cluster/projects/pughlab/projects/M4/ConsensusCruncher/processed/all_bams/$1 \
    -o $output/$1
