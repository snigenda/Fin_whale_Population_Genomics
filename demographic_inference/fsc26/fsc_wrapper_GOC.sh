#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=3G,highp
#$ -pe shared 8
#$ -N fscWrapperGOC
#$ -o /u/project/rwayne/pnunez/reports/fscwrapperGOC.out.txt
#$ -e /u/project/rwayne/pnunez/reports/fscwrapperGOC.err.txt
#$ -M pnunez
#$ -t 1-100

# This is a wrapper that will run 100 fastsimcoal iterations for each population for any list of models
# Author: Anabell Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub fsc_wrapper.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail


# Defined directories -------------

wd=/u/home/p/pnunez/project-rwayne/FinWhale/fastsimcoal
infDir=$wd/fastsimcoal_inference # specify where inference is happening
genericDir=$wd/ModelsFiles # location of generic FSC models
sfsDir=/u/home/p/pnunez/project-rwayne/FinWhale/SFS/SFS_projection_fsc
Â# sfsDir=/u/home/p/pnunez/project-rwayne/FinWhale/fastsimcoal/optimization_tryout

# Programs -------------------------

fsc=/u/home/p/pnunez/programs/fsc26_linux64/fsc26

# Parameters -----------------------

models='1D.4Epoch.GOC'
pops="GOC" 
cores=8 
ss="30"
version="version_normal_3"

########################################  MAIN  #############################################


for pop in $pops
do

	# Get sample size for the population 

	for model in $models
	do
		
		echo "starting $pop, $model"
		header=${model}

		# Copy generic files into directory and update 
		outdir=$infDir/$pop/$model/$version/run_${SGE_TASK_ID} 
		mkdir -p $outdir 

		cp $genericDir/$model.tpl $genericDir/$model.est $outdir # copy .est and .tpl files to outdir
		sed -i'' "s/SAMPLE_SIZE/$ss/g" $outdir/$model.tpl # sub in the sample size; note you need double quotes for variable to be expanded
		
		# Get sfs into inference directory and rename to match .est and .tpl files 
		cp $sfsDir/${pop}_MAFpop0.obs $outdir/${header}_MAFpop0.obs # copy your sfs into the directory where you'll be doing the fsc inference 
		cd $outdir
		$fsc -t ${header}.tpl -n 1000000 -m -e ${header}.est -M -L 60 -c${cores} -q


	done

done



