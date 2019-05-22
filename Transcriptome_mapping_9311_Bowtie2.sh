#!/bin/bash -e

# Bash script to map 9311 (minus & plus NAI) libraries to newly generated
# 9311 transcriptome reference
# For generating 9311 index from 9311 transcriptome refer script:
# index_builder_9311_ref_transcriptome.sh

#SBATCH --job-name=9311_transcriptome_mapping_bwt2
#SBATCH -o Batch_Bowtie2_9311_map.out
#SBATCH -e Batch_Bowtie2_9311_map.err
#SBATCH --mem 64000
#SBATCH -p nbi-long   # Partition (queue equivalent)


source bowtie2-2.2.9
source samtools-1.5


# map -NAI library, biologial replicate 1
bowtie2 --no-unal --no-mixed -p 4 -k 8 -x ./Index/9311trans -1 ./9mNAI_Bio1_forwardread_trimmed.fastq.gz -2 ./9mNAI_Bio1_reverseread_trimmed.fastq.gz -S 9mNAI_B1.sam
samtools view -bhS 9mNAI_B1.sam > 9mNAI_B1.bam
samtools sort -n 9mNAI_B1.bam -o 9mNAI_B1.srtname.bam 


# map -NAI library, biologial replicate 2
bowtie2 --no-unal --no-mixed -p 4 -k 8 -x ./Index/9311trans -1 ./9mNAI_Bio2_forwardread_trimmed.fastq.gz -2 ./9mNAI_Bio2_reverseread_trimmed.fastq.gz -S 9mNAI_B2.sam
samtools view -bhS 9mNAI_B2.sam > 9mNAI_B2.bam
samtools sort -n 9mNAI_B2.bam -o 9mNAI_B2.srtname.bam 


# map +NAI library, biologial replicate 1
bowtie2 --no-unal --no-mixed -p 4 -k 8 -x ./Index/9311trans -1 ./9pNAI_Bio1_forwardread_trimmed.fastq.gz -2 ./9pNAI_Bio1_reverseread_trimmed.fastq.gz -S 9pNAI_B1.sam 
samtools view -bhS 9pNAI_B1.sam > 9pNAI_B1.bam
samtools sort -n 9pNAI_B1.bam -o 9pNAI_B1.srtname.bam 


# map +NAI library, biologial replicate 2
bowtie2 --no-unal --no-mixed -p 4 -k 8 -x ./Index/9311trans -1 ./9pNAI_Bio2_forwardread_trimmed.fastq.gz -2 ./9pNAI_Bio2_reverseread_trimmed.fastq.gz -S 9pNAI_B2.sam 
samtools view -bhS 9pNAI_B2.sam > 9pNAI_B2.bam
samtools sort -n 9pNAI_B2.bam -o 9pNAI_B2.srtname.bam
