#
SNP calling pipeline readme and troubleshooting
Data for 93-11 (Indica) Minus libraries for SHAPE
Date: 12th January 2018

##################################################


#################################################
# Step 5. SNP Analysis
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ Prerequisites:									@@
@	jre-1.8.*	(using version 1.9 crashes GATK)				@@
@	samtools-1.5									@@
@	picard-1.134									@@
@	Latest version GATK3.8 jar file [if problems use nightly build]			@@
@					Also see Troubleshooting section		@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/SNP_analysis/

      5A. Generate reference .dict and .fai file for input Nippponbare genomic fasta reference
	# /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/References 

	(i) To get .dict file, run command

        Make sure no spaces in the fasta header (identifier). Spaces causes problems in samtools faidx

	source jre-1.8.0_45
	nohup java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar CreateSequenceDictionary R=Nipponbare_ref_genome.fasta O=Nipponbare_ref_genome.dict &

        NOTE:   The input fasta reference file has to have extension .fasta (Nipponbare_ref_genome.fasta not Nipponbare_ref_genome.con, else Picard complains!)
		The output extension has to be the input reference fasta file name and extension .dict (example should be 'Nipponbare_ref_genome.dict' not 'Nipponbare_ref_genome.fasta.dict', else GATK complains!)

        (ii) To get .fai file, run command

	source samtools-1.5
        samtools faidx ./Nipponbare_ref_genome.fasta

        The .fai output file (Nipponbare_ref_genome.fasta.fai) will be placed in the same directory where the input reference .fasta file is present.


				================================================
				| Run the pipeline			   	|
				| sbatch GATK_SNPcalling_pipeline.sh 		|
				|					   	|
				| OR individually run the following scripts: 	|
				|						|
				| NOTE: When running pipeline crosscheck if 	|
				| markduplicates (step 5B.v) crashes when 	|
				| running for the first time.			|
				|						|
				================================================


     5B. Use HaplotypeCaller from GATK
        # /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/SNP_analysis/

        Since this protocol is for RNAseq mapped to genome we use the guidelines mentioned here in steps 2 to 7:
	https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq


	(iii) Add Read Group (@RG) headers
              First, check if there are ReadGroup headers in the .bam file. If not, add @RG headers.
	      Sorting is optional here if the input BAM file is already sorted by position/co-ordinate.
	      Run script 'sbatch SLURM_addRG.sh'

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



	(iv) Merge 2 BAM files (Biological replicates B1 and B2)
	     Run script 'sbatch SLURM_mergeBAMfiles.sh' 

		srun samtools merge 9mNAI_B1B2.addRG.merged.bam 9mNAI_B1.addRG.bam 9mNAI_B2.addRG.bam



	 (v) Mark Duplicates
	     (@@ see Troubleshooting @@)
             Next, mark duplicate reads (skipping optical duplicates).
	     Run script 'sbatch SLURM_markDuplicates.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar MarkDuplicates \
			I=9mNAI_B1B2.addRG.merged.bam \
			O=9mNAI_B1B2.addRG.merged.dedupped.bam \
			READ_NAME_REGEX=null \
			CREATE_INDEX=true \
			VALIDATION_STRINGENCY=SILENT \
			M=output.metrics



	(vi) Split'N'Trim
	     Trim reads by hardclipping any regions overhanging into the intronic regions.
	     Run script 'sbatch SLURM_splitNTrim.sh' (@@ see Troubleeshooting @@)

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T SplitNCigarReads \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-I 9mNAI_B1B2.addRG.merged.dedupped.bam \
			-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.bam \
			-U ALLOW_N_CIGAR_READS 



       (vii) Indel Realignment
	     Local realignment around indels to correct mapping errors made by genomic aligners
	     (@@ see Troubleshooting @@)
	     Involves 2 steps:

	     (a) Make target intervals
		 Run script 'sbatch SLUM_indelRealigner-TargetCreator.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T RealignerTargetCreator \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.bam \
			-known ./../Indelrealignment/Nipponbare_indel.vcf \
			-nt 4 \
			-o 9mNAI_B1B2.realignertargetcreator.intervals


	     (b) Realign indels
		 Run script 'sbatch SLURM_indelRealigner.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.bam \
			-known ./../Indelrealignment/Nipponbare_indel.vcf \
			-targetIntervals 9mNAI_B1B2.realignertargetcreator.intervals \
		 	-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.bam


      (viii) BaseQualityScoreRecalibration (BQSR)
	     (@@ see Troubleshooting @@)
	     Involves 2 steps:

	     (a) Make recalibration table
		 Run script 'sbatch SLURM_BQSR_recaltable.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.bam \
			-knownSites /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/References/NB_bialSNP_pseudo_canonical_ALL.vcf \
			-o 9mNAI_B1B2_BaseRecalibrationReport_vcf_data.table
	    

	     (b) Recalibrate BAM file
		 Run script 'sbatch SLURM_BQSR.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T PrintReads \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.bam \
			-BQSR 9mNAI_B1B2_BaseRecalibrationReport_vcf_data.table \
			-nct 4 \
			-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.bam


	
	(ix) Single Nucleotide Variant (SNV) calling to get first set of RAW VARIANTS (unfiltered SNPs)
	     Run script 'sbatch SLURM_HaploTypeCaller.sh' (@@ see Troubleshooting @@)

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-I 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.bam \
			-maxAltAlleles 10 \
			-dontUseSoftClippedBases \
			-stand_call_conf 20.0 \
			-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-VARIANTS.vcf



	 (x) Variant filtering
	     Select high-confidence SNPs from the RAW VARIANTS using hard thresholds on specific annotations
	     Involves 2 steps:

	     (a) Extract only SNPs from the raw variant file
		 Run script 'sbatch SLURM_selectVariants.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T SelectVariants \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-V 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-VARIANTS.vcf \
			-selectType SNP \
			-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-SNPs.vcf 

	     (b) Apply hard filters to mark good quality SNPs
		 Run script 'sbatch SLURM_variantFilter.sh'

		srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/GenomeAnalysisTK-nightly-2018-01-01-1/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R ./../References/Nipponbare_ref_genome.fasta \
			-V 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.RAW-SNPs.vcf \
			--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
			--filterName "my_snp_filter" \
			-o 9mNAI_B1B2.addRG.merged.dedupped.splitntrimmed.indelrealigned.BQSR.FILTERED-SNPs.vcf
	 
		

###################################################################









@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        TROUBLE SHOOTING
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

NOTE	Use Samtools-1.5
	Use jre-1.8 (java)


1.      ** ERROR MESSAGE: SAM/BAM/CRAM file /nbi/Research-Groups/JIC/Yiliang-Ding/YD/ANALYSIS/OSA/HONG-RICE/Trial-SNP/9mNAI-B1/minus_alignment.bam is malformed.
	Please see https://software.broadinstitute.org/gatk/documentation/article?id=1317for more information. Error details: SAM file doesn't have any read groups defined in the header.  The GATK no longer supports SAM files without read groups**

        Therefore add @RG header lines and sort the bam file by coordinates:
        https://gatkforums.broadinstitute.org/gatk/discussion/2909
        https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        https://software.broadinstitute.org/gatk/documentation/article.php?id=59

        TO do this: run the SLURM_addRG.sh script
        srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar AddOrReplaceReadGroups INPUT=input.bam OUTPUT=output.bam \
		RGID=FCHKMVJBBXX_L2 \
		RGLB=9mNAI_B1 \
        	RGPL=illumina \
        	RGPU=FCHKMVJBBXX.L2.CGATGT \
		RGSM=9mNAI_B1 \
		SORT_ORDER=coordinate\
		CREATE_INDEX=true


2.      While running GATK with the safe version GATK3.8 causes java IO ERROR
        ** ERROR StatusLogger Log4j2 could not find a logging implementation. Please add log4j-core to the classpath. Using SimpleLogger to log to the console...
        **

        To fix this use the jar file from the nightly builds 2018-01-01-1

        https://gatkforums.broadinstitute.org/gatk/discussion/10004/realignertargetcreator-hangs
        https://software.broadinstitute.org/gatk/download/nightly


3.	To check the mapping quality use Picard CollectAlignmentSummaryMetrics
	srun java -jar /nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/Picard/picard.jar CollectAlignmentSummaryMetrics \
		INPUT=input.bam \
		OUTPUT=Picard_AlignmentSummaryMetric.out \
		R=./../References/Nipponbare_ref_genome.fasta
	
        The percentage of aligned reads for paired is SUM(PF_HQ_ALIGNED_READS)/PF_READS


4.	Running SLURM_markDuplicates.sh for the first time sometimes crashes due to a java IO error
	If this is the case re-run again and this time it works fine

	** ERROR Exception in thread "main" htsjdk.samtools.SAMException: /tmp/CSPI.8498237925943242181.tmp/28738.tmpnot found
        	at htsjdk.samtools.util.FileAppendStreamLRUCache$Functor.makeValue(FileAppendStreamLRUCache.java:64)
	        at htsjdk.samtools.util.FileAppendStreamLRUCache$Functor.makeValue(FileAppendStreamLRUCache.java:49)
	        at htsjdk.samtools.util.ResourceLimitedMap.get(ResourceLimitedMap.java:76)
	        at htsjdk.samtools.CoordinateSortedPairInfoMap.getOutputStreamForSequence(CoordinateSortedPairInfoMap.java:180)
	        at htsjdk.samtools.CoordinateSortedPairInfoMap.put(CoordinateSortedPairInfoMap.java:164)
	        at picard.sam.markduplicates.util.DiskBasedReadEndsForMarkDuplicatesMap.put(DiskBasedReadEndsForMarkDuplicatesMap.java:65)
	        at picard.sam.markduplicates.MarkDuplicates.buildSortedReadEndLists(MarkDuplicates.java:535)
	        at picard.sam.markduplicates.MarkDuplicates.doWork(MarkDuplicates.java:232)
	        at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:268)
	        at picard.cmdline.PicardCommandLine.instanceMain(PicardCommandLine.java:98)
	        at picard.cmdline.PicardCommandLine.main(PicardCommandLine.java:108)
	Caused by: java.io.FileNotFoundException: /tmp/CSPI.8498237925943242181.tmp/28738.tmp (Too many open files)
	        at java.io.FileOutputStream.open0(Native Method)
	        at java.io.FileOutputStream.open(FileOutputStream.java:270)
	        at java.io.FileOutputStream.<init>(FileOutputStream.java:213)
	        at htsjdk.samtools.util.FileAppendStreamLRUCache$Functor.makeValue(FileAppendStreamLRUCache.java:61)
	        ... 10 more
	srun: error: n128n21: task 0: Exited with exit code 1


5.	Filtering optical duplicates while marking duplicates (Step 5B.v) has no effect on the number of reads filtered by HCMappingQualityFilter
	but affects Library Size Estimation.


6.	To generate an intervals VCF file to use for Indel Realignment (step 5B.vii.a)
	(By Hongjing)
	Download BED, BIM and FAM files of 3K Rice Genomes 
	http://snp-seek.irri.org/_download.zul under "3K RG 2.3mio biallelic indel Dataset"/"Subsets 3k RG SNPs release 1.0".
	To convert it to .vcf following https://www.biostars.org/p/108499/ and https://www.cog-genomics.org/plink/1.9/data#recode
	The command line is:
	plink --bfile Nipponbare_indel --recode vcf --out Nipponbare_indel

	The resulting .vcf file has only chromosome number in 1st coloumn. Modify the first coloumn to include Chr1 for example.
	This is required because the output SAM and BAM files have Chr1 and not 1 or Chr01.

	sed -i 's/^/Chr/g' Nipponbare_indel.vcf
	sed -i 's/^Chr#/#/g' Nipponbare_indel.vcf


7.	To use known goldstandard vcf files for generating base recalibration table (step 5B.viii.a)
	Download VCF file from http://snp-seek.irri.org/_download.zul under "Effect of 29mio biallelic SNPs on Rice Genome Annotation Project rel 7 gene models"/"SnpEff results for 3K RG 29mio biallelic SNPs, VCF file"
	
	The resulting .vcf file has only chromosome number in 1st coloumn. Modify the first coloumn to include Chr1 for example.
	This is required because the output SAM and BAM files have Chr1 and not 1 or Chr01.
	
	sed -i 's/^/Chr/g' NB_bialSNP_pseudo_canonical_ALL.
	sed -i 's/^Chr#/#/g' NB_bialSNP_pseudo_canonical_ALL.vcf
