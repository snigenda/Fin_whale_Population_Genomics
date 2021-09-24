#!/bin/bash
#$ -l h_data=12G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/ROH/combine_bcftools_roh_20210104.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/ROH/combine_bcftools_roh_20210104.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub step0_combine_bcftools_20210104.sh
# @description	Combine the bcftools output to a shortened version
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Jan  6 17:27:11 2021

###########################################################
## import packages

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
WORKSCRIPT=${HOMEDIR}/scripts/important_results/Runs_of_homozygosity/step0_combine_bcftools_20210104.R

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

Rscript --vanilla ${WORKSCRIPT}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"