#!/bin/bash -e

## Bash script to submit job for calculating nucleotide conservation scores

#SBATCH --job-name=SNPcons
#SBATCH -o SNPcons.out
#SBATCH -e SNPcons.err
#SBATCH --mem 500gb
#SBATCH-p nbi-largemem   # Partition (queue equivalent)


## Usage: sbatch SLURM_get_SNPconservation_3Kgenomes.sh 


source perl-5.22.1

perl get_snp_conservation_3Kgenomes.pl Mapped_Rice3K_accessions_varieties_types.out 3kall_snpposition_map.tsv Universe_matrix_geno_NB
