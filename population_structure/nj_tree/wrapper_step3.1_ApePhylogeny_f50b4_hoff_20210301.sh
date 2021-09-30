#!/bin/bash
#$ -l highp,h_data=12G,h_vmem=16G,h_rt=72:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step3.1_ApePhylogeny_f50b4_hoff_20210301.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step3.1_ApePhylogeny_f50b4_hoff_20210301.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub wrapper_step3_ApePhylogeny_f50b4_20210301.sh <gdsfile> <outprefix>
# @description	Wrapper to submit the Rscript for Ape Trees with bootstraps (NOTE: Using three maf cutoffs) Not plotting anything
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Mar  2 01:57:00 2021

############################################################
## import packages

sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
GDSFILE=${1}
OUTPREFIX=${2}
DATASET='f50b4'
REF='Minke'
TODAY=$(date "+%Y%m%d")

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}
LOGDIR=${WORKDIR}/logs
mkdir -p ${LOGDIR}

WORKSCRIPT=${HOMEDIR}/scripts/PopStructure/f50b4/step3.1_ApePhylogeny_f50b4_hoff_20210301.R
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

LOG=${LOGDIR}/step3.1_ApePhylogeny_${OUTPREFIX}_${TODAY}.log

############################################################
## main

cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; GDSFILE = ${GDSFILE}; OUTPREFIX = ${OUTPREFIX}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; GDSFILE = ${GDSFILE}; OUTPREFIX = ${OUTPREFIX}" > ${LOG}

Rscript --vanilla ${WORKSCRIPT} ${WORKDIR} ${GDSFILE} ${OUTPREFIX} &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}

conda deactivate