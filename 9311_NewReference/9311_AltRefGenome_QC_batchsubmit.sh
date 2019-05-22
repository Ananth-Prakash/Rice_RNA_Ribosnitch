#!/bin/bash

# Bash script to run all vcf files in parallel for QC of newly generated 93-11 Genome

RES_1=/nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/Alternate_93-11_Reference
RES_2=/nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Hongjing_Rice/HISAT2_Mapping_Genomic/93-11_Genomic_mapping/References
CWD=$( pwd )

cd ${CWD}

mapfile -t INPUTS < input_vcf_file_names.list

sleep 10

#echo ${INPUTS[${SLURM_PROCID}]}
echo ./${INPUTS[${SLURM_PROCID}]}.vcf

source perl-5.22.1

perl Genome_slice_QC.pl ${CWD}/${INPUTS[${SLURM_PROCID}]}.vcf ${RES_1}/93-11_AltReferenceGenome.fasta ${RES_2}/Nipponbare_ref_genome.fasta > ${CWD}/ReferenceFasta_GenomeQC_${INPUTS[${SLURM_PROCID}]}.out

