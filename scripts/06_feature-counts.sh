#!/bin/bash

#SBATCH -p batch
#SBATCH -J counts
#SBATCH -n 8

#Read count quantification
#Script to count the number of reads aligned to the reference genome using featureCounts.

module load featureCounts

cd /home/emurungi/gitau/marion/raw

# process paired-end data

BAM_DIR="/home/emurungi/gitau/marion/TNBC/TNBCbam/"

bam_file_name=$(basename "$BAM_DIR" .sorted.bam)
GTF_FILE=GCF_000001405.39_GRCh38.p13_genomic.gtf


for bam_file in ${BAM_DIR}/*.sorted.bam; do

	featureCounts -p -O -T 8 -a $GTF_FILE -o ${bam_file}.counts.txt $bam_file

done