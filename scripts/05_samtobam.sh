#!/bin/bash

#SBATCH -p batch
#SBATCH -J sam-to-bam
#SBATCH -n 8

# load the blast module
module load samtools

mkdir -p /home/gitau/marion/TNBC/TNBCbam

cd /home/gitau/marion/TNBC/TNBCsam

SAM_DIR="/home/gitau/marion/TNBC/TNBCsam"

BAM_DIR="/home/gitau/marion/TNBC/TNBCbam/"


# convert sam file to sorted bam files
for sam_file in ${SAM_DIR}/*.sam; do
        sam_file_name=$(basename "$sam_file" .sam)
        samtools view -S -b $sam_file | \
        samtools sort -n -o ${BAM_DIR}/${sam_file_name}.sorted.bam
done

