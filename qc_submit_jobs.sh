#!/bin/bash -e

## Script to assess quality of the reads [QC]

#SBATCH --job-name=QC
#SBATCH -o Batch_QC.out
#SBATCH -e Batch_QC.err
#SBATCH --mem 2000
#SBATCH -p nbi-short   # Partition (queue equivalent)


source fastqc-0.11.3

# Nipponbare
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B1_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B1_R2.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B2_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NmNAI_B2_R2.fastq.gz -o ./'

sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B1_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B1_R2.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B2_R1.fastq.gz -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/NpNAI_B2_R2.fastq.gz -o ./'


# 9311
sbatch --wrap 'fastqc ./../Trimmed_files/9mNAI_Bio1_forwardread_trimmed.fastq -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/9mNAI_Bio1_reverseread_trimmed.fastq -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/9mNAI_Bio2_forwardread_trimmed.fastq -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/9mNAI_Bio2_reverseread_trimmed.fastq -o ./'

sbatch --wrap 'fastqc ./../Trimmed_files/9pNAI_Bio1_forwardread_trimmed.fastq -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/9pNAI_Bio1_reverseread_trimmed.fastq -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/9pNAI_Bio2_forwardread_trimmed.fastq -o ./'
sbatch --wrap 'fastqc ./../Trimmed_files/9pNAI_Bio2_reverseread_trimmed.fastq -o ./'
