#!/bin/bash -e

# Step 5B.iii in README.txt. 
# Before Variant Calling using GATK. (https://gatkforums.broadinstitute.org/gatk/discussion/2909)
#				     (https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq)
# 				     (https://software.broadinstitute.org/gatk/documentation/article.php?id=59)
# To add missing Read Group headers "@RG"


#SBATCH --job-name=Picard_addRG
#SBATCH -o Picard_addRG_9mNAI.out
#SBATCH -e Picard_addRG_9mNAI.err
#SBATCH --mem 64GB
#SBATCH -n 1
#SBATCH -p nbi-long   # Partition (queue equivalent)

source jre-1.8.0_45 

# /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/Genomic_mapping_HISAT2

# For 9mNAI-B1.
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar AddOrReplaceReadGroups \
	INPUT=./../Genomic_mapping_HISAT2/9mNAI_B1.srt.bam \
	OUTPUT=9mNAI_B1.addRG.bam \
	RGID=FCHKMVJBBXX_L2 \
	RGLB=9mNAI_B1 \
	RGPL=illumina \
	RGPU=FCHKMVJBBXX.L2.CGATGT \
	RGSM=9mNAI_B1 \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true


# For 9mNAI-B2
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar AddOrReplaceReadGroups \
	INPUT=./../Genomic_mapping_HISAT2/9mNAI_B2.srt.bam \
	OUTPUT=9mNAI_B2.addRG.bam \
	RGID=FCHKMVJBBXX_L3 \
	RGLB=9mNAI_B2 \
	RGPL=illumina \
	RGPU=FCHKMVJBBXX.L3.TTAGGC \
	RGSM=9mNAI_B2 \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true
