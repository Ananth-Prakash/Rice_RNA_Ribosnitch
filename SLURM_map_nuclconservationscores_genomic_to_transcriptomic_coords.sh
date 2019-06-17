#!/bin/bash -e

# Bash script to genomic positions in the output file of 3K conservation (got from running SLURM_get_SNPconservation_3Kgenomes.sh)
# scores into transcript positions 

#SBATCH --job-name=SNPposition_map_genomic_to_transc
#SBATCH -o SNPposition_map.out
#SBATCH -e SNPposition_map.err
#SBATCH --mem 32GB
#SBATCH -n 1
#SBATCH -p nbi-long   # Partition (queue equivalent)


source perl-5.22.1


perl map_nuclconservationscores_genome_transcriptome_coords.pl all_gff3_MSU7.gff3 SNPConservation_All_3K_Accessions.out > SNPConservation_All_3K_Accessions_transcript_positions_mapped.out
