#!/bin/bash -e

## Shell script to map NAI minus libraries of 9311
## biological replicates to Nipponbare reference
## genome.
 

#SBATCH --job-name=hisat2_9311_minus
#SBATCH -o Batch_HISAT2_9mNAI_map.out
#SBATCH -e Batch_HISAT2_9mNAI_map.err
#SBATCH --mem 64GB
#SBATCH -p nbi-long   # Partition (queue equivalent)


source switch-institute ei
source HISAT-2.0.5
source samtools-1.5


## 9311 Minus NAI library Biological replicate 1

srun hisat2 --rna-strandness FR -k 8 --no-unal -p 4 -x ./../Genomic_mapping_HISAT2/Index/Nipgenome -1 ./9mNAI_Bio1_forwardread_trimmed.fastq.gz -2 ./9mNAI_Bio1_reverseread_trimmed.fastq.gz -S 9mNAI_B1.sam 
srun samtools view -bhS 9mNAI_B1.sam > 9mNAI_B1.bam
srun samtools sort 9mNAI_B1.bam -o 9mNAI_B1.srt.bam



## 9311 Minus NAI library Biological replicate 2

srun hisat2 --rna-strandness FR -k 8 --no-unal -p 4 -x ./../Genomic_mapping_HISAT2/Index/Nipgenome -1 ./9mNAI_Bio2_forwardread_trimmed.fastq.gz -2 ./9mNAI_Bio2_reverseread_trimmed.fastq.gz -S 9mNAI_B2.sam 
srun samtools view -bhS 9mNAI_B2.sam > 9mNAI_B2.bam
srun samtools sort 9mNAI_B2.bam -o 9mNAI_B2.srt.bam
