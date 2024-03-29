#!/bin/bash

#SBATCH -p batch
#SBATCH -J QC
#SBATCH -n 8

module load fastqc

#mkdir -p /home/gitau/marion/results/TNBCfastqc

#Script to generate FastQC reports using the FastQC tool

# fastq files directory
FASTQ_DIR="/home/gitau/marion/TNBC/TNBCtrimmed"

# fastqc reports directory
REPORT_DIR="/home/gitau/marion/TNBC/fastqc-report"

for file in $FASTQ_DIR/*.fastq; do
   fastqc ${file} -o ${REPORT_DIR}
done
