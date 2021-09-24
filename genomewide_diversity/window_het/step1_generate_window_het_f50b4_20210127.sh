#!/bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/window_het/step1_generate_window_het_f50b4_20210127.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/window_het/step1_generate_window_het_f50b4_20210127.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub -t 1-<ncontig> step1_generate_window_het_f50b4_20210127.sh
# @description	Generate windowed heterozygosity along contigs, allows multiple contigs in one vcf file (for the baleen_genomes `f50b4` dataset)
# Author: Jacqueline Robinson, Sergio Nigenda, Meixi Lin (meixilin@ucla.edu)
# Date: Wed Jan 27 13:46:50 2021

###########################################################
## import packages
sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
DATASET="f50b4"
REF="Minke"
TODAY=$(date "+%Y%m%d")
IDX=$(printf %02d ${SGE_TASK_ID})

# directories
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
OUTDIR=${HOMEDIR}/window_het/${DATASET}/${REF}/window_het_${TODAY}
LOGDIR=${OUTDIR}/logs
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

# workscript and log
WORKSCRIPT=${HOMEDIR}/scripts/window_het/SlidingWindowHet_finwhale.py
LOG=${LOGDIR}/step1_generate_window_het_${IDX}_${TODAY}.log

# vcf info
VCFDIR=${HOMEDIR}/baleen_genomes/filteredvcf/${DATASET}/${REF}
VCFPREFIX="JointCalls_${DATASET}_08_B_VariantFiltration_${IDX}"

# window het settings
WINDOWSIZE=1000000
STEPSIZE=1000000

###########################################################
## main
cd ${OUTDIR}
# logging
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; git commit id: ${COMMITID}" > ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] COMMAND: python ${WORKSCRIPT} \
--path "${OUTDIR}" \
--VCF "${VCFDIR}/${VCFPREFIX}.vcf.gz" \
--outprefix "${VCFPREFIX}" \
--windowsize "${WINDOWSIZE}" \
--stepsize "${STEPSIZE}" \
--idx "${IDX}"" >> ${LOG}

# start analyses
python ${WORKSCRIPT} \
--path "${OUTDIR}" \
--VCF "${VCFDIR}/${VCFPREFIX}.vcf.gz" \
--outprefix "${VCFPREFIX}" \
--windowsize "${WINDOWSIZE}" \
--stepsize "${STEPSIZE}" \
--idx "${IDX}" &>> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}

conda deactivate
