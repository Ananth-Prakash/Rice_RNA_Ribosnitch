#!/bin/bash -e

# Bash script to run trimmomatic pipleine

#SBATCH --job-name=maticM
#SBATCH -o Batch_maticM.out
#SBATCH -e Batch_maticM.err
#SBATCH --mem 2000
#SBATCH -p nbi-short   # Partition (queue equivalent)

export SBATCH_PARTITION=nbi-long 

source  trimmomatic-0.33 
#source cutadapt-1.9.1 

#This script works for NmNAI (-) files
srun job_trim1_Nipon.sh 
srun job_trim2_Nipon.sh 
