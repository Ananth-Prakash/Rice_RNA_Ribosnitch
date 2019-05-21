#!/bin/bash -e

## Pipeline to call SNPs in rice population using Nipponbare genome as reference
## The pipeline follows GATK best practices for SNP calling
## Troubleshooting at each step is documented in Readme. 

#SBATCH --job-name=pipeline_GATK_SNP
#SBATCH -o pipeline_GATK_SNP.out
#SBATCH -e pipeline_GATK_SNP.err
#SBATCH --mem 64GB
#SBATCH -n 1
#SBATCH -p nbi-long   # Partition (queue equivalent)

source jre-1.8.0_45


echo -e "["$(date)"]\tMarking duplicates.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar MarkDuplicates \
	I=9mNAI_B1B2.addRG.merged.bam \
	O=9mNAI_B1B2.addRG.merged.dedupped.bam \
	READ_NAME_REGEX=null \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT \
	M=output.metrics 2>MarkDuplicates.log

echo -e "["$(date)"]\tSpliting reads.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T SplitNCigarReads \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-I 9mNAI_B1B2.addRG.merged.dedupped.bam \
	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.bam \
	-U ALLOW_N_CIGAR_READS 2>SplitNCigarReads.log

echo -e "["$(date)"]\tCreating targets for indel realignment.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.bam \
	-known ./../Indelrealignment/Nipponbare_indel.vcf \
	-nt 4 \
 	-o 9mNAI_B1B2.realignertargetcreator.intervals 2>indeltarget.log

echo -e "["$(date)"]\tPerforming Indel Realignment.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.bam \
	-known ./../Indelrealignment/Nipponbare_indel.vcf \
	-targetIntervals 9mNAI_B1B2.realignertargetcreator.intervals \
 	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.bam 2>indelrealigned.log

echo -e "["$(date)"]\tPerforming BQSR.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.bam \
	-knownSites /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/References/NB_bialSNP_pseudo_canonical_ALL.vcf \
	-o 9mNAI_B1B2_BaseRecalibrationReport_vcf_data.table 2>BQSRtable.log

echo -e "["$(date)"]\tPrinting recalibrated reads.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T PrintReads \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.bam \
	-BQSR 9mNAI_B1B2_BaseRecalibrationReport_vcf_data.table \
	-nct 4 \
	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.bam 2>BQSRprint.log

echo -e "["$(date)"]\tRunning HaplotypeCaller.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.bam \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-VARIANTS.vcf 2>HaplotypeCaller.log

echo -e "["$(date)"]\tSelecting SNPs.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-V 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-VARIANTS.vcf \
	-selectType SNP \
	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-SNPs.vcf 2>SelectSNPs.log

echo -e "["$(date)"]\tFiltering Variants.."
java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R ./../References/Nipponbare_ref_genome.fasta \
	-V 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-SNPs.vcf \
	--filterExpression "QD < 2.0 || FS > 30.0" \
	--filterName "my_snp_filter" \
	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.FILTERED-SNPs.vcf 2>VariantFilter.log
