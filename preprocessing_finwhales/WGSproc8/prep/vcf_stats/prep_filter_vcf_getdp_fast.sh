#! /bin/bash
#$ -l h_data=12G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -o /u/project/rwayne/snigenda/finwhale/reports/prep_filter_vcf_getdp_all50.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/prep_filter_vcf_getdp_all50.err.txt
#$ -m abe

# Usage: generate per contiglist sum on gtDP for all sites
# Please specify the amount of iteration needed based on the DNA_FILE genome
# qsub -t 1-23 prep_filter_vcf_stats_fast.sh meixilin Bryde
# qsub -t 1-96 prep_filter_vcf_stats_fast.sh meixilin Minke
# Sun May 17 20:27:11 2020
# @modification Tue Aug 11 13:43:07 2020
# @modification adapt for all50 dataset
# @modification: Fri Nov  6 20:42:09 2020
# @modification: Rerun to match the SIFT update

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
# INPUT
USER=${1}
REF=${2}
# Variables related to contig list
IDX=$(printf %02d ${SGE_TASK_ID})

# Directories
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKDIR=${HOMEDIR}/filteredvcf/all50/${REF}
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale
OUTDIR=${WORKDIR}/variant_summary
LOGDIR=${WORKDIR}/logs
LOG=${LOGDIR}/0x_${IDX}_${REF}_get_gtdp.log

# Variables
TODAY=$(date "+%Y%m%d")
COMMITID=`git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master`
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc0/utils/extract_vcf_gtdp.py

###########################################################
## main

# logging
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${SGE_TASK_ID} Generate statistics file all50 ${REF} ${IDX} using ${WORKSCRIPT}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${SGE_TASK_ID}" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Generate statistics file all50 ${REF} ${IDX} using ${WORKSCRIPT}" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${LOG}

mkdir -p ${OUTDIR}
cd ${OUTDIR}

python ${WORKSCRIPT} \
--VCF ${WORKDIR}/JointCalls_all50_07_snpEff_SIFT_formatted_${IDX}.vcf.gz \
--outfile "${OUTDIR}/all50_${REF}_chr${IDX}_${TODAY}" \
--contiglist "$IDX" &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${SGE_TASK_ID} Done generating statistics"
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${SGE_TASK_ID} Done generating statistics" >> ${LOG}

conda deactivate