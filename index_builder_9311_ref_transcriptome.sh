#!/bin/bash -e

# Bash script to run index builder of 9311 transcriptome fasta
# for mapping 9311 reads to 9311 transcriptome

#SBATCH --job-name=index
#SBATCH -o Batch_Bowtie2_index_build.out
#SBATCH -e Batch_Bowtie2_index_build.err
#SBATCH --mem 32000
#SBATCH -p nbi-short   # Partition (queue equivalent)


source bowtie2-2.3.1

bowtie2-build -f 9311.transcriptome.fasta 9311trans
