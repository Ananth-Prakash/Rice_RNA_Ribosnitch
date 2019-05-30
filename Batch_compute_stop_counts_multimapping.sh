#!/bin/bash -e

# Bash script to run counting script in batch mode

#SBATCH --job-name=counts
#SBATCH -o Batch_stopcounts.out
#SBATCH -e Batch_stopcounts.err
#SBATCH --mem 64000
#SBATCH -p nbi-medium   # Partition (queue equivalent)

source python-3.5.1

# Nipponbare
python Compute_stop_countsM.py NmNAI_B1.srtname.bam NmNAI_B1.cnt
python Compute_stop_countsM.py NmNAI_B2.srtname.bam NmNAI_B2.cnt
python Compute_stop_countsM.py NpNAI_B1.srtname.bam NpNAI_B1.cnt
python Compute_stop_countsM.py NpNAI_B2.srtname.bam NpNAI_B2.cnt

# 9311
python Compute_stop_countsM.py 9mNAI_B1.srtname.bam 9mNAI_B1.cnt
python Compute_stop_countsM.py 9mNAI_B2.srtname.bam 9mNAI_B2.cnt
python Compute_stop_countsM.py 9pNAI_B1.srtname.bam 9pNAI_B1.cnt
python Compute_stop_countsM.py 9pNAI_B2.srtname.bam 9pNAI_B2.cnt

