#!/bin/bash
#$ -l h_data=24G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step2_Prunedgds2vcf_all50_20210316.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step2_Prunedgds2vcf_all50_20210316.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub wrapper_step2_Prunedgds2vcf_all50_20210316.sh
# @description	Wrapper to submit the Rscript for converting LD pruned gds to vcf files (NOTE: Using three maf cutoffs)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Mar 18 14:16:53 2021

############################################################
## import packages

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
DATASET='all50'
REF='Minke'
MAFLISTS=('05' '10' 'NA')

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}

WORKSCRIPT=${HOMEDIR}/scripts/PopStructure/all50/step2_Prunedgds2vcf_all50_20210316.R
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}"

for MAFCUT in ${MAFLISTS[@]}; do
echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Start."
# The input files
OUTPREFIX=JointCalls_${DATASET}_filterpass_bialleic_all_LDPruned_maf${MAFCUT}
GDSFILE=${OUTPREFIX}.gds
Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX}
echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Done."
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate