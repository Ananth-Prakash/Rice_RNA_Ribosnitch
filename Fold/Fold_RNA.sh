#!/bin/bash -e

# Script to fold RNA sequences:
# (i)  SHAPE constrained (invivo folding) by providing input SHAPE data
# (ii) SHAPE unconstrained (insilico folding) without any inputSHAPE data
# Reads file names to fold from input list.

#SBATCH --job-name=Fold_RNA
#SBATCH -o Fold_RNA.out
#SBATCH -e Fold_RNA.err
#SBATCH --mem 8gb
#SBATCH-p nbi-medium   # Partition (queue equivalent)


## Usage: sbatch Fold_RNA.sh file.list.inp


source viennaRNA-2.3.3

# assign input & output file suffixes and prefixes
fasta_file='.fasta'
shape_file='.SHAPE'
shape_prefix='m1.9b-0.5_'
unconstr_prefix='UNCONSTR_'

# output directory for SHAPE constrained folds
mkdir SHAPE_CONSTR_folds
# output directory for SHAPE unconstrained folds
mkdir UNCONSTR_folds

while IFS='' read -r line || [[ -n "$line" ]]
do

	inp_fasta=$line$fasta_file
	inp_shape=$line$shape_file
	shape_constr_out_fold=$shape_prefix$line
	shape_unconstr_out_fold=$unconstr_prefix$line

	# (i)  SHAPE constrained folding (invivo)
	RNAfold -p -d2 -T 28.0 -i FASTA_files/$inp_fasta  --shapeMethod="Dm1.9b-0.5" --shape=SHAPE_files/$inp_shape  --outfile=SHAPE_CONSTR_folds/$shape_constr_out_fold

	# (ii) SHAPE unconstrained folding (insilico) 
	RNAfold -p -d2 -T 28.0 -i FASTA_files/$inp_fasta  --outfile=UNCONSTR_folds/$shape_unconstr_out_fold
	
done < "$1"
