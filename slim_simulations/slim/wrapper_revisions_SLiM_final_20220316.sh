#!/bin/bash
#$ -l highp,h_rt=140:00:00,h_data=20G,h_vmem=20G,arch=intel-gold*
#$ -wd <homedir>
#$ -o <homedir>/reports/revisions_SLiM_final/wrapper_revisions_SLiM_final_20220316.out.txt
#$ -e <homedir>/reports/revisions_SLiM_final/wrapper_revisions_SLiM_final_20220316.err.txt
#$ -m abe
#$ -t 1-25

# @version      v1
# @usage        qsub wrapper_revisions_SLiM_final_20220316.sh <prefix> <full path to slim worker file>
# @description  wrapper to submit any given slim file during revisions. Final SLiM simulations to be included in the main text.
# Author: Meixi Lin
# Date: 2022-03-16 11:52:34

###########################################################
## import packages
sleep $((RANDOM % 30))

set -o pipefail

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

###########################################################
## def functions

###########################################################
## def variables
TODAY=$(date "+%Y%m%d")
PREFIX=${1}
JOBFILE=${2}

SLIM=/u/project/klohmuel/ckyriazi/software/slim_build/slim # software to use
HOMEDIR=<homedir>
WORKDIR=${HOMEDIR}/revisions_SLiM_final/${PREFIX}/rep${SGE_TASK_ID}
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

OUTFILE=${WORKDIR}/slimfinal_${PREFIX}_rep${SGE_TASK_ID}.txt
LOG=${WORKDIR}/slimfinal_${PREFIX}_rep${SGE_TASK_ID}.log

if [ -d ${WORKDIR} ]; then
    rm -r ${WORKDIR}
fi

###########################################################
## main
mkdir -p ${WORKDIR}
cd ${WORKDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}" > ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] COMMAND: $SLIM ${JOBFILE} > ${OUTFILE}"
echo -e "[$(date "+%Y-%m-%d %T")] COMMAND: $SLIM ${JOBFILE} > ${OUTFILE}" >> ${LOG}

$SLIM ${JOBFILE} > ${OUTFILE}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}

