#!/bin/bash

# Script to process fold file of a transcript to calculate
# basepairing information.

#SBATCH --job-name=Relplot
#SBATCH -o Relplot_log.out
#SBATCH -e Relplot_log.err
#SBATCH --mem 8gb
#SBATCH-p nbi-long   # Partition (queue equivalent)


## Usage: sbatch Relplot.sh RNAfold_output_files_ss.ps.list

# Before running Replot fold RNA using RNAfold.
# RNAfold_output_files_ss.ps.list file consists of list of 
# UNIQUE FILE NAMES WITHOUT FILE EXTENSIONS that are the
# output of RNAfold alogrithm
# The RNAfold output files for a transcript are 
# 	*_ss.ps, 
#	*_dp.ps and 
#	*.fold 


# Add input and output files suffix and prefix
ss_ext='_ss.ps'
dp_ext='_dp.ps'
out_ext='_rss.ps'


while IFS='' read -r line || [[ -n "$line" ]]
do
	inp_ss=$line$ss_ext
	inp_dp=$line$dp_ext

	out_rel=$line$out_ext

	/nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/ViennaRNA-2.4.6/src/Utils/relplot.pl -a ./RNAfold_output/$inp_ss ./RNAfold_output/$inp_dp > ./Relplot_output/$out_rel
	
done < $1

