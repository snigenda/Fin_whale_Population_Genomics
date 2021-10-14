#! /bin/bash
#$ -cwd
#$ -l h_rt=48:00:00,h_data=3G,highp
#$ -pe shared 8
#$ -N fscWrapper
#$ -o /u/project/rwayne/pnunez/reports/fscwrapper.out.txt
#$ -e /u/project/rwayne/pnunez/reports/fscwrapper.err.txt
#$ -M pnunez
#$ -t 1-100

# This is a wrapper that will run 100 fastsimcoal iterations for joint SFS for any list of models
# Author: Anabell Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub fsc_wrapper.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail


# Defined directories -------------

wd=/u/home/p/pnunez/project-rwayne/FinWhale/fastsimcoal
infDir=$wd/fastsimcoal_inference # specify where inference is happening
genericDir=$wd/ModelsFiles # location of generic FSC models
sfs=/u/home/p/pnunez/project-rwayne/FinWhale/SFS/SFS_projection_fsc/ENP_GOC_jointMAFpop1_0.obs

# Programs -------------------------

fsc=/u/home/p/pnunez/programs/fsc26_linux64/fsc26

# Parameters -----------------------

models='2D.2Epoch 2D.2Epoch.mig 2D.3Epoch 2D.3Epoch.mig'
cores=8
version="version_3"

########################################  MAIN  #############################################


for model in $models
do

        echo "starting $model"
        header=${model}

        # Copy generic files into directory and update
        outdir=$infDir/2DModels/$model/$version/run_${SGE_TASK_ID}
        mkdir -p $outdir

        cp $genericDir/$model.tpl $genericDir/$model.est $outdir # copy .est and .tpl files to outdir

        # Get sfs into inference directory and rename to match .est and .tpl files
        cp $sfs $outdir/${header}_jointMAFpop1_0.obs # copy your sfs into the directory where you'll be doing the fsc inference
        cd $outdir
        $fsc -t ${header}.tpl -n 1000000 -m -e ${header}.est -M -L 60 -c${cores} -q


done

