#!/bin/bash -e

#SBATCH --job-name=maticM
#SBATCH -o Batch_maticM.out
#SBATCH -e Batch_maticM.err
#SBATCH --mem 2000
#SBATCH -p nbi-short   # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-long 

source fastqc-0.11.3

sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B1_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B1_R2.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B2_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B2_R2.fastq.gz -o ./'

sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B1_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B1_R2.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B2_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B2_R2.fastq.gz -o ./'
