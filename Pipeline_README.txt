
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
		[python generate-TranscriptReference-usingAssemblySeq.py]
------------


## STEP 7. Map 9311 reads (+NAI and -NAI) to new 9311 Transcriptome
		Use Bowtie2

	   (i)	Create index for transcriptomic mapping with 9311.transcriptome.fasta
		Run script
		[sbatch index_builder_9311_ref_transcriptome.sh]

	  (ii)	Map (-)NAI and (+)NAI libraries to new reference of 93-11
		Run script
		[sbatch Transcriptome_mapping_9311_Bowtie2.sh]
------------


## STEP 8. Computing stop counts and coverage

	   (i)	Run script
		[sbatch Batch_compute_stop_counts_multimapping.sh]
		output is a .cnt file which has (i) Transcript id, (ii) stopcounts and (ii) coverage at each nucleotide position.

	  (ii)	Summarise the stopcounts and coverage and normalise by transcript length
		using the .cnt file obtained from running the first step
		Run script
		[perl stopcounts_coverage_fraction.pl input.cnt > output.cvg]
------------


# STEP 9. Calculate SHAPE reactivity

	   (i)	Merge (sum) the counts from two biological replicates.
        	activate rnaenv
        	rna_structure counts -l 9mNAI_B1.cnt 9mNAI_B2.cnt 9mNAI_merged.cnt
        	rna_structure counts -l 9pNAI_B1.cnt 9pNAI_B2.cnt 9pNAI_merged.cnt

	  (ii)	Calculate the box plot normalized reactivity
		rna_structure_cli shape-reactivity -v --calc=log-noalpha --norm=boxplot-noq3 --norm-exclude-zeroes --out-filter=pos --cut-nts=-40 9311_combinedREF.fasta 9mNAI_merged.cnt 9pNAI_merged.cnt ./Boxplot_react/ trans > transcriptome_bpnorm.log 2>&1

	 (iii)	Calculate 95% Winsorization reactivity of passed files
		Run script
        	[perl Rescale_reactivity_Winsorize.pl Raw_reactivity_pass.list Rescaled_react_pass_directory/]
------------


# STEP 10. Fold RNA secondary structures
	
	   (i)  Fold Nipponbare an 9311 transcripts
		SHAPE constrained folding (invivo folding) and unconstrained folding (insilico)
		Provide input list of file names without extensions (fasta shape files must have the same names) to fold
		Run script
		[sbatch Fold_RNA.sh file.list.inp]
	
	 (ii) 	Annotate basepairing partners
		Identify pairing partners later to calculate basepairing distances and probabilities 
		Provide list of name_ss.ps files from RNAfold as input
		Run script
		[sbatch Relplot.sh RNAfold_output_files_ss.ps.list]
	       
	(iii)	Extract basepairing probabilities
		Provide input files, which are obtained from running Relplot.
		Run script
		[perl Calc_BasePairProbability.pl path_to_fold_rss.ps_files.list]
------------


# STEP 11. Compare RNA secondary structures

	  (i)	PPV


------------


# STEP 12. Calculate nucleotide conservation score

		The multiallelic data of 32 Million SNP positions in 3000 (3K) rice accessions was downloaded from
		https://s3.amazonaws.com/3kricegenome/reduced/3k_RG_32mio_All_multiallelic_biallelic_SNP_dataset.zip
		OR the header ("3K RG 32mio SNPs, called vs Nipponbare MSU7/IRGSP1.0 genome, tabular format") in http://snp-seek.irri.org/_download.zul

	  (i)	Prepare mapping of rice IRIS_UNIQUE_IDS to Variety Group (ex. Indica, Japonica, etc.)
		Download rice Variety Group information of 3K rice accessions from
		"The 3,000 rice genomes project"
		Additional file 1: Table S1A:
		https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4035669/#S1
		
		Extract columns C ("DNA_UNIQUE_ID") and O ("Variety Group (Tree)1") in sheet "Table-1A_IRRI" &
			columns C ("DNA_UNIQUE_ID") and N ("Variety Group (Tree)2") in sheet "Table-1B_CAAS"
		concatenate the coloumns into a file.

	 	Run script
		[Map_3KRice_accession_variety_name.R]

	 (ii)	Compute conservation scores of nucleotides
		Run script
		[sbatch SLURM_get_SNPconservation_3Kgenomes.sh]

	(iii)	Map nucleotide conservation score from genomic cordinates to transcript cordinates
		Run script
		[sbatch SLURM_map_nuclconservationscores_genomic_to_transcriptomic_coords.sh]
