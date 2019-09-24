#!/bin/bash

#Script to extract/compute basepairing partners from fold file 

#SBATCH --job-name=Relplot_SHAPE
#SBATCH -o Relplot_SHAPE.out
#SBATCH -e Relplot_SHAPE.err
#SBATCH --mem 8gb
#SBATCH-p nbi-medium   # Partition (queue equivalent)


## Usage: sbatch Relplot.sh

export SBATCH_PARTITION=nbi-medium

source viennaRNA-2.3.3

#list of ss files
# SHAPE_CONSTR_FullLength/*ss.ps
# UNCONSTR_FullLength/*ss.ps 


sed -i -e 's/_ss.ps$//g' shape_constr_inp_files


ss_ext='_ss.ps'
dp_ext='_dp.ps'
out_ext='_rss.ps'


while IFS='' read -r line || [[ -n "$line" ]]
do
	inp_ss=$line$ss_ext
	inp_dp=$line$dp_ext

	out=$(echo $line | sed -e 's/.*\///g')
	out_rel=$out$out_ext

	/nbi/Research-Groups/JIC/Yiliang-Ding/ANANTH/Software_p/ViennaRNA-2.4.6/src/Utils/relplot.pl -a $inp_ss $inp_dp > ./Relplot/SHAPE_CONSTR_FullLength/$out_rel
	
done < shape_constr_inp_files

