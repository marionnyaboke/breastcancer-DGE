#!/bin/bash

#SBATCH -p batch
#SBATCH -J multiqc
#SBATCH -n 8

module load multiqc

#Script to generate MultiQC reports using the MultiQC tool

mkdir -p /home/emurungi/gitau/marion/results/multiqc

# fastqc directory
FASTQC_DIR="/home/emurungi/gitau/marion/results/fastqc1"

# Multiqc output directory
OUT_DIR="/home/emurungi/gitau/marion/results/multiqc"

# run mutiqc
multiqc ${FASTQC_DIR}/*_fastqc.zip --outdir ${OUT_DIR}
