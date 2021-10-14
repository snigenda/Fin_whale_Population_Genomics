#! /bin/bash
#$ -cwd
#$ -l h_rt=00:09:00,h_data=8G
#$ -N gatherResults
#$ -o /u/project/rwayne/pnunez/reports/fscresults.out.txt
#$ -e /u/project/rwayne/pnunez/reports/fscresults.err.txt
#$ -m abe
#$ -M pnunez

# Gather up FSC results
# Author: Anabell Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub fsc_gather_results.sh


#Define directories ------------

wd=/u/home/p/pnunez/project-rwayne/FinWhale/fastsimcoal
infDir=$wd/fastsimcoal_inference # specify where inference is happening


#Define Parameters -----------

models='1D.3Epoch.ENP'
pops="ENP"
version="version_normal_10"


############################ MAIN ##############################

mkdir $infDir/resultsSummaries

for pop in $pops
do

	for model in $models
	do
		sumdir=$infDir/resultsSummaries/$pop/$model/$version                                                                               
		mkdir -p $sumdir

		header=${model}

		outfile=$sumdir/${model}.results.Sumaries.csv

		header=`head -n1 $infDir/$pop/$model/$version/run_1/$model/*bestlhoods`
		echo -e "runNum\t$header" | tr "\\t" "," > $outfile

		for i in {1..100}
		do
			outdir=$infDir/$pop/$model/$version/run_${i}/$model 
			results=`grep -v [A-Z] $outdir/*bestlhoods`
			echo -e "$i\t$results"| tr "\\t" "," >> $outfile

		done
		
		
	done

done
