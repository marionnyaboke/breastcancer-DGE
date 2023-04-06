#!/bin/bash

#SBATCH -p batch
#SBATCH -J alignment
#SBATCH -n 8

# load the module
module load hisat2

mkdir -p /home/gitau/marion/TNBC/TNBCsam

cd /home/gitau/marion/TNBC/TNBCtrimmed

FNA_DIR="/home/gitau/marion/raw"

SAM_DIR="/home/gitau/marion/TNBC/TNBCsam"

#Align fastqs to indexed reference genome
SAMPLES="SRR10729843 SRR10729844 SRR10729846 SRR10729847 SRR10729848 SRR10729849 SRR10729850 SRR10729851 SRR10729852 SRR10729853 SRR10729854 SRR10729855 SRR10729856 SRR10729857 SRR10729858"

for SAMPLE in $SAMPLES; do

        hisat2 \
                 -x ${FNA_DIR}/GCF_000001405.39_GRCh38.p13_genomic.fna_index_hisat2 \
                 -1 ${SAMPLE}_1.trimmed.fastq \
                 -2 ${SAMPLE}_2.trimmed.fastq \
                 -S ${SAM_DIR}/${SAMPLE}.sam \
                 -p 6 \
                --summary-file ${SAMPLE}.txt \
                --new-summary

done



