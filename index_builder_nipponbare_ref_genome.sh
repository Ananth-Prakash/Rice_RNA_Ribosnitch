#!/bin/bash -e

# Bash script to run index builder pipleine

#SBATCH --job-name=index
#SBATCH -o Batch_HISAT2_index_build.out
#SBATCH -e Batch_HISAT2_index_build.err
#SBATCH --mem 64000
#SBATCH -p nbi-medium   # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-long

source switch-institute ei
source HISAT-2.0.5

srun hisat2-build /nbi/Research-Groups/JIC/Yiliang-Ding/YD/ANALYSIS/OSA/HONG-RICE/Mapping-NAI-Nipp-20171113/all.chrs.con Nipgenome --ss ./../all.ss --exon ./../all.exon
