#!/bin/bash

#Trim script for Nm (-) files Bio 1

read1=./../Merged_160721_160918_1_NmNAI_1.fq.gz
read2=./../Merged_160721_160918_2_NmNAI_1.fq.gz

baseout=NmNAI_B1

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${read1} ${read2}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21


### this will the final output prefix
#baseout=mSHAPE-rep1

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_fw_paired.fq   ${baseout}_R1.fastq  HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_rev_paired.fq  ${baseout}_R2.fastq  HEADCROP:0




#Trim script for Nm (-) files Bio 2

read1=./../Merged_160721_160918_1_NmNAI_2.fq.gz
read2=./../Merged_160721_160918_2_NmNAI_2.fq.gz

baseout=NmNAI_B2

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${read1} ${read2}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21


### this will the final output prefix
#baseout=mSHAPE-rep1

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_fw_paired.fq   ${baseout}_R1.fastq  HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8  ${baseout}_out_rev_paired.fq  ${baseout}_R2.fastq  HEADCROP:0




#Trim script for Np (+) files Bio 1

read1=./../Merged_160721_160918_1_NpNAI_1.fq.gz
read2=./../Merged_160721_160918_2_NpNAI_1.fq.gz

baseout=NpNAI_B1

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${read1} ${read2}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21

# TRAILING:20 means low-quality nucleotides at the 3'-end (<20) will be trimmed.
## MINLEN:21 means any read with a length < 21 will not be included.


### this will the final output prefix
#baseout=pSHAPE-K

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8     ${baseout}_out_fw_paired.fq    ${baseout}_R1.fastq HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8     ${baseout}_out_rev_paired.fq   ${baseout}_R2.fastq HEADCROP:0 




#Trim script for Np (+) files Bio 2

read1=./../Merged_160721_160918_1_NpNAI_2.fq.gz
read2=./../Merged_160721_160918_2_NpNAI_2.fq.gz

baseout=NpNAI_B2

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar  PE -threads 8  ${read1} ${read2}  \
${baseout}_out_fw_paired.fq ${baseout}_out_fw_unpaired.fq ${baseout}_out_rev_paired.fq ${baseout}_out_rev_unpaired.fq \
ILLUMINACLIP:/nbi/software/testing/trimmomatic/0.33/x86_64/bin/adapters/TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:20  MINLEN:21

# TRAILING:20 means low-quality nucleotides at the 3'-end (<20) will be trimmed.
## MINLEN:21 means any read with a length < 21 will not be included.


### this will the final output prefix
#baseout=pSHAPE-K

java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8     ${baseout}_out_fw_paired.fq    ${baseout}_R1.fastq HEADCROP:3 
java -jar /nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar   SE -threads 8     ${baseout}_out_rev_paired.fq   ${baseout}_R2.fastq HEADCROP:0 



