#!/bin/bash
#$ -l h_data=24G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step1_gdsLDPruning_all50_20210316.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step1_gdsLDPruning_all50_20210316.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub wrapper_step1_gdsLDPruning_all50_20210316.sh
# @description	Wrapper to submit the Rscript for LDPruning (NOTE: Using three maf cutoffs)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Mar 17 11:26:56 2021

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

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}

WORKSCRIPT=${HOMEDIR}/scripts/PopStructure/all50/step1_gdsPASSLDPruning_all50_20210316.R
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# The input files
GDSFILE=JointCalls_all50_08_B_VariantFiltration_bialleic_all.gds
OUTPREFIX=JointCalls_all50_filterpass_bialleic_all

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}"

Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate