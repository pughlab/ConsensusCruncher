#/bin/bash
#$ -S /bin/bash
#$ -cwd

module load python3/3.4.3
module load bwa/0.7.2
module load samtools/0.1.18

python3 /mnt/work1/users/pughlab/src/ConsensusCruncher/consensus_scripts/extract_barcodes.py --read1 /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/fastq//LargeMid_59_L006_R1.fastq --read2 /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/fastq//LargeMid_59_L006_R2.fastq --outfile /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/fastq_tag/LargeMid_59_L006 --blen 2 --slen 1 --sfilt T 

bwa mem -M -t4 -R '@RG	ID:1	SM:LargeMid_59_L006	PL:Illumina	PU:59.L006	LB:ConsensusCruncher_tester' /mnt/work1/data/genomes/human/hg19/iGenomes/Sequence/BWAIndex/genome.fa /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/fastq_tag/LargeMid_59_L006'_barcode_R1.fastq' /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/fastq_tag/LargeMid_59_L006'_barcode_R2.fastq' > /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/bamfiles/LargeMid_59_L006.sam 

samtools view -bhS /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/bamfiles/LargeMid_59_L006.sam | samtools sort -@4 - /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/bamfiles/LargeMid_59_L006 

samtools index /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/bamfiles/LargeMid_59_L006.bam
rm /mnt/work1/users/pughlab/src/ConsensusCruncher/tester/bamfiles/LargeMid_59_L006.sam
