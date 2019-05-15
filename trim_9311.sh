#!/bin/bash

##
## Script to trim adapters of 9311 minus and plus libraries
##


#SBATCH --job-name=trim_9311
#SBATCH -o trim_9311.out
#SBATCH -e trim_9311.err
#SBATCH --mem 64GB
#SBATCH -p nbi-medium   # Partition (queue equivalent)

source  trimmomatic-0.33


##
## Trim adapters of 9311 minus NAI Biological Replicate 1

Forwardread=./../HD9mNAI-Bio1/Merged_forwardread_techreps_Bio1.fq.gz
Reverseread=./../HD9mNAI-Bio1/Merged_reverseread_techreps_Bio1.fq.gz

### this will add the final output prefix
baseout=9mNAI_Bio1

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${Forwardread} ${Reverseread}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21


java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_fw_paired.fq   ${baseout}_forwardread_trimmed.fastq  HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_rev_paired.fq  ${baseout}_reverseread_trimmed.fastq  HEADCROP:0



##
## Trim adapters of 9311 minus NAI Biological Replicate 2

Forwardread=./../HD9mNAI-Bio2/Merged_forwardread_techreps_Bio2.fq.gz
Reverseread=./../HD9mNAI-Bio2/Merged_reverseread_techreps_Bio2.fq.gz

### this will add the final output prefix
baseout=9mNAI_Bio2

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${Forwardread} ${Reverseread}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21


java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_fw_paired.fq   ${baseout}_forwardread_trimmed.fastq  HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_rev_paired.fq  ${baseout}_reverseread_trimmed.fastq  HEADCROP:0




##
## Trim adapters of 9311 plus NAI Biological Replicate 1

Forwardread=./../HD9pNAI-Bio1/Merged_forwardread_plus_Bio1.fq.gz
Reverseread=./../HD9pNAI-Bio1/Merged_reverseread_plus_Bio1.fq.gz

### this will add the final output prefix
baseout=9pNAI_Bio1

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${Forwardread} ${Reverseread}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21


java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_fw_paired.fq   ${baseout}_forwardread_trimmed.fastq  HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_rev_paired.fq  ${baseout}_reverseread_trimmed.fastq  HEADCROP:0



##
## Trim adapters of 9311 plus NAI Biological Replicate 2

Forwardread=./../HD9pNAI-Bio2/Merged_forwardread_plus_Bio2.fq.gz
Reverseread=./../HD9pNAI-Bio2/Merged_reverseread_plus_Bio2.fq.gz

### this will add the final output prefix
baseout=9pNAI_Bio2

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${Forwardread} ${Reverseread}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21


java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_fw_paired.fq   ${baseout}_forwardread_trimmed.fastq  HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_rev_paired.fq  ${baseout}_reverseread_trimmed.fastq  HEADCROP:0
