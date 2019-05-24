
## STEP 1.	Merge/Append [cat] technical replicates

# Oryza sativa sp. Japonica var. Nipponbare

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070025_1_NmNAI_1.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070025_1_NmNAI_1.fq.gz	(Forward strand technical replicates)	-NAI	Bio1
  

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070025_2_NmNAI_1.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070025_2_NmNAI_1.fq.gz	(Reverse strand technical replicates)	-NAI	Bio1
  ------------

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070027_1_NpNAI_1.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070027_1_NpNAI_1.fq.gz	(Forward strand technical replicates)	+NAI	Bio1

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070027_2_NpNAI_1.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070027_2_NpNAI_1.fq.gz	(Reverse strand technical replicates)	+NAI	Bio1
  ------------

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070026_1_NmNAI_2.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070026_1_NmNAI_2.fq.gz	(Forward strand technical replicates)	-NAI	Bio2

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070026_2_NmNAI_2.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070026_2_NmNAI_2.fq.gz	(Reverse strand technical replicates)	-NAI	Bio2
  ------------

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070028_1_NpNAI_2.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070028_1_NpNAI_2.fq.gz	(Forward strand technical replicates)	+NAI	Bio2

  160721_I114_FCHC2FYBBXX_L2_CHKPE85216070028_2_NpNAI_2.fq.gz
  160918_I211_FCHCGLCBBXX_L6_CHKPE85216070028_2_NpNAI_2.fq.gz	(Reverse strand technical replicates)	+NAI	Bio2


  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070025_1_NmNAI_1.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070025_1_NmNAI_1.fq.gz > Merged_160721_160918_1_NmNAI_1.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070025_2_NmNAI_1.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070025_2_NmNAI_1.fq.gz > Merged_160721_160918_2_NmNAI_1.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070027_1_NpNAI_1.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070027_1_NpNAI_1.fq.gz > Merged_160721_160918_1_NpNAI_1.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070027_2_NpNAI_1.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070027_2_NpNAI_1.fq.gz > Merged_160721_160918_2_NpNAI_1.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070026_1_NmNAI_2.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070026_1_NmNAI_2.fq.gz > Merged_160721_160918_1_NmNAI_2.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070026_2_NmNAI_2.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070026_2_NmNAI_2.fq.gz > Merged_160721_160918_2_NmNAI_2.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070028_1_NpNAI_2.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070028_1_NpNAI_2.fq.gz > Merged_160721_160918_1_NpNAI_2.fq.gz]
  [cat 160721_I114_FCHC2FYBBXX_L2_CHKPE85216070028_2_NpNAI_2.fq.gz 160918_I211_FCHCGLCBBXX_L6_CHKPE85216070028_2_NpNAI_2.fq.gz > Merged_160721_160918_2_NpNAI_2.fq.gz]

  Merged_160721_160918_1_NmNAI_1.fq.gz				(Merged & compressed forward strand replicates; Nipponbare -NAI Bio1)
  Merged_160721_160918_2_NmNAI_1.fq.gz				(Merged & compressed reverse strand replicates; Nipponbare -NAI Bio1)
  ------------
  Merged_160721_160918_1_NpNAI_1.fq.gz				(Merged & compressed forward strand replicates; Nipponbare +NAI Bio1)
  Merged_160721_160918_2_NpNAI_1.fq.gz				(Merged & compressed reverse strand replicates; Nipponbare +NAI Bio1)
  ------------
  Merged_160721_160918_1_NmNAI_2.fq.gz				(Merged & compressed forward strand replicates; Nipponbare -NAI Bio2)
  Merged_160721_160918_2_NmNAI_2.fq.gz				(Merged & compressed reverse strand replicates; Nipponbare -NAI Bio2)
  ------------
  Merged_160721_160918_1_NpNAI_2.fq.gz				(Merged & compressed forward strand replicates; Nipponbare +NAI Bio2)
  Merged_160721_160918_2_NpNAI_2.fq.gz				(Merged & compressed reverse strand replicates; Nipponbare +NAI Bio2)
  ------------


# Oryza sativa sp. Indica var. 9311

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090001_1.fq.gz		
  	       FCHKMVJBBXX_L2_CHKPE85216090001_1.fq.gz		(Forward strand technical replicates) 	-NAI	Bio1

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090001_2.fq.gz		
 	       FCHKMVJBBXX_L2_CHKPE85216090001_2.fq.gz		(Reverse strand technical replicates) 	-NAI	Bio1
  -----------------

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090003_1.fq.gz		
	       FCHKMVJBBXX_L2_CHKPE85216090003_1.fq.gz		(Forward strand technical replicates)	+NAI	Bio1

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090003_2.fq.gz		
	       FCHKMVJBBXX_L2_CHKPE85216090003_2.fq.gz		(Reverse strand technical replicates)	+NAI	Bio1
  -----------------

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090002_1.fq.gz		
	       FCHKMVJBBXX_L3_CHKPE85216090002_1.fq.gz		(Forward strand technical replicates)	-NAI	Bio2

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090002_2.fq.gz		
  	       FCHKMVJBBXX_L3_CHKPE85216090002_2.fq.gz		(Reverse strand technical replicates)	-NAI	Bio2
  -----------------

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090004_1.fq.gz		
	       FCHKMVJBBXX_L3_CHKPE85216090004_1.fq.gz		(Forward strand technical replicates)	+NAI	Bio2

  161003_I211_FCHCH53BBXX_L1_CHKPE85216090004_2.fq.gz		
	       FCHKMVJBBXX_L3_CHKPE85216090004_2.fq.gz		(Reverse strand technical replicates)	+NAI	Bio2
  -----------------

  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090001_1.fq.gz FCHKMVJBBXX_L2_CHKPE85216090001_1.fq.gz > Merged_forwardread_techreps_Bio1.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090001_2.fq.gz FCHKMVJBBXX_L2_CHKPE85216090001_2.fq.gz > Merged_reverseread_techreps_Bio1.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090003_1.fq.gz FCHKMVJBBXX_L2_CHKPE85216090003_1.fq.gz > Merged_forwardread_plus_Bio1.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090003_2.fq.gz FCHKMVJBBXX_L2_CHKPE85216090003_2.fq.gz > Merged_reverseread_plus_Bio1.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090002_1.fq.gz FCHKMVJBBXX_L3_CHKPE85216090002_1.fq.gz > Merged_forwardread_techreps_Bio2.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090002_2.fq.gz FCHKMVJBBXX_L3_CHKPE85216090002_2.fq.gz > Merged_reverseread_techreps_Bio2.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090004_1.fq.gz FCHKMVJBBXX_L3_CHKPE85216090004_1.fq.gz > Merged_forwardread_plus_Bio2.fq.gz]
  [cat 161003_I211_FCHCH53BBXX_L1_CHKPE85216090004_2.fq.gz FCHKMVJBBXX_L3_CHKPE85216090004_2.fq.gz > Merged_reverseread_plus_Bio2.fq.gz]

  Merged_forwardread_techreps_Bio1.fq.gz			(Merged & compressed forward strand replicates; 9311 -NAI Bio1)
  Merged_reverseread_techreps_Bio1.fq.gz			(Merged & compressed reverse strand replicates; 9311 -NAI Bio1)
  ------------
  Merged_forwardread_plus_Bio1.fq.gz				(Merged & compressed forward strand replicates; 9311 +NAI Bio1)
  Merged_reverseread_plus_Bio1.fq.gz				(Merged & compressed forward strand replicates; 9311 +NAI Bio1)
  ------------
  Merged_forwardread_techreps_Bio2.fq.gz			(Merged & compressed forward strand replicates; 9311 -NAI Bio2)
  Merged_reverseread_techreps_Bio2.fq.gz			(Merged & compressed reverse strand replicates; 9311 -NAI Bio2)
  ------------
  Merged_forwardread_plus_Bio2.fq.gz				(Merged & compressed forward strand replicates; 9311 +NAI Bio2)
  Merged_reverseread_plus_Bio2.fq.gz				(Merged & compressed forward strand replicates; 9311 +NAI Bio2)
  ------------


## STEP 2. Adapter trimming using Trimmomatic

# Nipponbare
  Run script 
   [sbatch trim_nipponbare.sh]

  NmNAI_B1_R1.fastq.gz		(NAI minus_BiologicalRep. 1_Forward read)
  NmNAI_B1_R2.fastq.gz		(NAI minus_BiologicalRep. 1_Reverse read)
  NmNAI_B2_R1.fastq.gz		(NAI minus_BiologicalRep. 2_Forward read)
  NmNAI_B2_R2.fastq.gz		(NAI minus_BiologicalRep. 2_Reverse read)
  NpNAI_B1_R1.fastq.gz		(NAI plus_BiologicalRep. 1_Forward read)
  NpNAI_B1_R2.fastq.gz		(NAI plus_BiologicalRep. 1_Reverse read)
  NpNAI_B2_R1.fastq.gz		(NAI plus_BiologicalRep. 2_Forward read)
  NpNAI_B2_R2.fastq.gz		(NAI plus_BiologicalRep. 2_Reverse read)

# 9311
  Run script
   [sbatch trim_nipponbare.sh]

  9mNAI_B1_forwardread_trimmed.fastq.gz	(NAI minus_BiologicalRep. 1_Forward read)
  9mNAI_B1_reverseread_trimmed.fastq.gz	(NAI minus_BiologicalRep. 1_Reverse read)
  9mNAI_B2_forwardread_trimmed.fastq.gz	(NAI minus_BiologicalRep. 2_Forward read)
  9mNAI_B2_reverseread_trimmed.fastq.gz	(NAI minus_BiologicalRep. 2_Reverse read)
  9pNAI_B1_forwardread_trimmed.fastq.gz	(NAI plus_BiologicalRep. 1_Forward read)
  9pNAI_B1_reverseread_trimmed.fastq.gz	(NAI plus_BiologicalRep. 1_Reverse read)
  9pNAI_B2_forwardread_trimmed.fastq.gz	(NAI plus_BiologicalRep. 2_Forward read)
  9pNAI_B2_reverseread_trimmed.fastq.gz	(NAI plus_BiologicalRep. 2_Reverse read)

The output files *_out_fw_paired.fq
	  	 *_out_fw_unpaired.fq
	  	 *_out_rev_paired.fq
	  	 *_out_rev_unpaired.fq		are deleted.

Output files	 *.fastq.gz
	  	 *.fastq.gz			are zipped and kept for QC and mapping (Step 3 and 4).
------------


## STEP 3. QC using FastQC
   Run script
   [sbatch qc_submit_jobs.sh]
------------


## STEP 4. Mapping minus libraries (-NAI) of 9311 to Nipponbare reference genome
  	   9311 has no reference genome (or transcriptome) to map reads of 9311 to, 
 	   therefore to generate a reference 9311 genome (through SNP calling on Nipponbare genome and replacement of SNPs), 
	   the minus 9311 libraries (9mNAI_B1 and 9mNAI_B2) are mapped to the Nipponbare genome.

   Run script
   [sbatch Genomic_mapping_9311_HISAT2.sh]
------------


## STEP 5. SNP calling
   Run script
   [sbatch GATK_SNPcalling_pipeline.sh]
------------


## STEP 6. Generate Reference Genome & Transcriptome for 9311

	   (i)	Generate reference Genome
 
	   	Use GATK FastaAlternateReferenceMaker  
		Run script
		[sbatch SLURM_AlternateReferenceMaker.sh]

	  (ii)	Quality check the newly generated 9311 genome by comparing it with Nipponbare reference genome
		Check if SNPs are incorporated in their correct positions in the new 9311 genome.
 
		Run script 
		[perl Genome_slice_QC.pl]
		or batch submit separately for each chromosome vcf
		[sbatch 9311_AltRefGenome_QC_batchsubmit.sh]

	 (iii)	Generate reference Transcriptome
		Run script
		[python generate-reference-usingAssemblySeq-IsoChk-final.py]
------------


## STEP 7. Map 9311 reads (+NAI and -NAI) to new 9311 Transcriptome
		Use Bowtie2

	   (i)	Create index for transcriptomic mapping with 9311.transcriptome.fasta
		Run script
		[sbatch index_builder_9311_ref_transcriptome.sh]

	  (ii)	Map (-)NAI and (+)NAI libraries to new reference of 93-11
		Run script
		[sbatch Transcriptome_mapping_9311_Bowtie2.sh]

