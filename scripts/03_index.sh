#!/bin/bash

#SBATCH -p batch
#SBATCH -J index
#SBATCH -n 8

# load the blast module
module load hisat2

cd /home/gitau/marion/raw

FNA_DIR="/home/gitau/marion/raw"

#Index reference using HISAT2
hisat2-build ${FNA_DIR}/GCF_000001405.39_GRCh38.p13_genomic.fna GCF_000001405.39_GRCh38.p13_genomic.fna_index_hisat2


