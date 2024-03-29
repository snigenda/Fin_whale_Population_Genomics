#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=23:00:00
#$ -wd <homedir2>/finwhale
#$ -o <homedir2>/finwhale/reports/WGSproc9_all50.out.txt
#$ -e <homedir2>/finwhale/reports/WGSproc9_all50.err.txt
#$ -m abe

# @version 		v0
# @usage		generate filter statistics
# @description	WGSproc9_all50
# Author: Meixi Lin
# Date: Fri Aug 14 10:58:02 2020
# @modification Mon Aug 24 22:00:34 2020
# @modification Final version

###########################################################
## import packages
sleep $((RANDOM % 120))


conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
DATASET="all50" # all48;testsix
REF="Minke" # Minke;Bryde
TODAY=$(date "+%Y%m%d")

HOMEDIR=<homedir2>/finwhale

WORKDIR=${HOMEDIR}/filteredvcf/${DATASET}/${REF}
SCRATCHDIR=<scratchdir>/finwhale/filteredvcf/${DATASET}/${REF}
STATDIR=<homedir>/Summary_stats/${DATASET}/${REF}/filter_stats_${TODAY}
mkdir -p ${SCRATCHDIR}
mkdir -p ${STATDIR}
mkdir -p ${STATDIR}/logs

IDX=$(printf %02d ${SGE_TASK_ID}) # this is SGE specific array id, essentially the contig list
MYPREFIX="${DATASET}_${REF}_${IDX}_filter_stats"
COMMITID=`git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master`

PARSESCRIPT="<homedir2>/finwhale/scripts/WGSproc9/parse_customVCFfilter_summary_20201110.py"
LOG=${STATDIR}/logs/${MYPREFIX}_${TODAY}.log

###########################################################
## main
cd ${SCRATCHDIR}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Filter stats"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Start generate filter stats for ${MYPREFIX}" >> ${LOG}

# parse and tally the filter statistics
python ${PARSESCRIPT} --dir "${SCRATCHDIR}" --prefix "${MYPREFIX}" --contig "${IDX}" --outdir "${STATDIR}"
echo -e "[$(date "+%Y-%m-%d %T")] ${PARSESCRIPT} Done" >> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}

conda deactivate
