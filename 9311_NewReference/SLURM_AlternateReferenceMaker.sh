#!/bin/bash -e

# Step 6 (in Pipeline)
# Generating Alternate Genomic Reference for 93-11 based on high confidence filtered SNPs
# https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_fasta_FastaAlternateReferenceMaker.php

#SBATCH --job-name=GATK_AltRef
#SBATCH -o GATK_AlternateReference.out
#SBATCH -e GATK_AlternateReference.err
#SBATCH --mem 32GB
#SBATCH -n 1
#SBATCH -p nbi-medium   # Partition (queue equivalent)

source jre-1.8.0_45

java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-V ./../SNP_analysis/9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.FILTERED-SNPs.vcf \
	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.FILTERED-SNPs-forref.vcf \
	-restrictAllelesTo BIALLELIC \
	--removeUnusedAlternates

java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T FastaAlternateReferenceMaker \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-V 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.FILTERED-SNPs-forref.vcf \
	-o 93-11_AltReferenceGenome.fasta
#	-L input.intervals \
#	[--snpmask mask.vcf]
