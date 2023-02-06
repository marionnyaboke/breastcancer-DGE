#!/bin/bash

#SBATCH -p batch
#SBATCH -J trimming
#SBATCH -n 8

module load cutadapt

#Script to trim adaptors


mkdir -p /home/emurungi/gitau/marion/TNBC/TNBCtrimmed


FASTQ_DIR="/home/emurungi/gitau/marion/TNBC/raw"

OUT_DIR="/home/emurungi/gitau/marion/TNBC/TNBCtrimmed"


SAMPLES="SRR10729843 SRR10729844 SRR10729846 SRR10729847 SRR10729848 SRR10729849 SRR10729850 SRR10729851 SRR10729852 SRR10729853 SRR10729854 SRR10729855 SRR10729856 SRR10729857 SRR10729858"

for SAMPLE in $SAMPLES; do

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m15 -o ${OUT_DIR}/${SAMPLE}_1.trimmed.fastq -p ${OUT_DIR}/${SAMPLE}_2.trimmed.fastq ${FASTQ_DIR}/${SAMPLE}_1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_2.fastq.gz

done


