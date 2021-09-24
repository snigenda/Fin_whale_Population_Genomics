#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc7_all50.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc7_all50.err.txt
#$ -m abe

# Step 7: Add annotations for the mutations using snpEff
# Adapted from Best Practices for GATK3
# Authors: Sergio & Meixi
# @modification Mon Aug 17 13:44:25 2020
# @modification Add SIFT option
# @modification Thu Nov  5 15:24:32 2020
# @modification 1. Rerun the WGSproc7 to update the formats and not backing up to sirius
# @modification 2. Add git commit ID to each run's report

# Examples: qsub -t 1-96/1-23 WGSproc7_a_snpEff_SIFT_annotations_20201105.sh meixilin Minke/Bryde

sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

################################################################################

### Set variables

USER=${1}
REF=${2}

BAMHEAD="MarkDuplicates"
SNPEFFDIR=/u/project/rwayne/software/finwhale/miniconda2/envs/gentools/share/snpeff-4.3.1t-2
SIFT=/u/project/rwayne/software/finwhale/sift/SIFT4G_Annotator.jar
HOMEDIR=/u/project/rwayne/snigenda/finwhale
# SIRIUSDIR=/data3/finwhale
COMMITID=`git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master`

WORKDIR=${HOMEDIR}/filteredvcf/all50/${REF}
# WORKSIRIDIR=${SIRIUSDIR}/filteredvcf/all50/${REF}
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale/filteredvcf/all50/${REF}

mkdir -p ${WORKDIR}
mkdir -p ${SCRATCHDIR}
# ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${WORKSIRIDIR}"

IDX=$(printf %02d ${SGE_TASK_ID}) # this is SGE specific array id, essentially the contig list
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc7/formatSIFTvcf_all50_20201105.py

# define variables by references
if [ $REF == 'Minke' ]; then
    SNPEFFDB=Baac01.10776
    SIFTDB=/u/project/rwayne/jarobins/utils/programs/sift/databases/GCF_000493695.1_BalAcu1.0
fi

if [ $REF == 'Bryde' ]; then
    SNPEFFDB=Baed01.141314
fi

################################################################################

### logging
# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Start WGSproc7 for all50 ${REF}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
cd ${WORKDIR}

# starts logging
LOGDIR=${WORKDIR}/logs/annotations
mkdir -p ${LOGDIR}
PROGRESSLOG=${WORKDIR}/logs/WGSproc7_${REF}_${BAMHEAD}_${IDX}_progress_all50.log
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] snpEff annotations based on ${SNPEFFDB} ... " >> ${PROGRESSLOG}

LOG=${LOGDIR}/04_A_all50_${REF}_${BAMHEAD}_snpEff_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}


################################################################################

### snpEFF annotations

snpEff -Xmx8g -nodownload -v -canon -stats ${LOGDIR}/${REF}_chr${IDX}.html \
${SNPEFFDB} \
JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz > \
${SCRATCHDIR}/JointCalls_all50_07_snpEff_${IDX}.vcf 2> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}
date "+%Y-%m-%d %T" >> ${LOG}

################################################################################

### SIFT annotations
LOG=${LOGDIR}/04_B_all50_${REF}_${BAMHEAD}_SIFT_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] SIFT annotations based on ${SIFTDB} ... " >> ${PROGRESSLOG}

mkdir -p ${SCRATCHDIR}/${IDX}
cd ${SCRATCHDIR}/${IDX}

java -Xmx8g -jar ${SIFT} \
-c -t -r ./ -d ${SIFTDB} -i ${SCRATCHDIR}/JointCalls_all50_07_snpEff_${IDX}.vcf &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
mv JointCalls_all50_07_snpEff_${IDX}_SIFTpredictions.vcf JointCalls_all50_07_snpEff_SIFT_${IDX}.vcf
rm ${SCRATCHDIR}/JointCalls_all50_07_snpEff_${IDX}.vcf

################################################################################

### format SIFT annotations
echo -e "[$(date "+%Y-%m-%d %T")] Formatting JointCalls_all50_07_snpEff_SIFT_${IDX}.vcf using ${WORKSCRIPT} ... " >> ${PROGRESSLOG}

LOG=${LOGDIR}/04_C_all50_${REF}_${BAMHEAD}_formatSIFT_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

python ${WORKSCRIPT} JointCalls_all50_07_snpEff_SIFT_${IDX}.vcf | \
bgzip > JointCalls_all50_07_snpEff_SIFT_formatted_${IDX}.vcf.gz
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

tabix -p vcf JointCalls_all50_07_snpEff_SIFT_formatted_${IDX}.vcf.gz
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}
date "+%Y-%m-%d %T" >> ${LOG}

# delete files
echo -e "[$(date "+%Y-%m-%d %T")] Cleaning up ... " >> ${PROGRESSLOG}
set -eo pipefail
rm JointCalls_all50_07_snpEff_SIFT_${IDX}.vcf

# move files and logs
mv Events.log ${LOGDIR}/Events_SIFTannotations_${IDX}.log
mv JointCalls_all50_07_snpEff_${IDX}_SIFTannotations.xls ${LOGDIR}/JointCalls_all50_07_snpEff_SIFT_${IDX}.xls
mv JointCalls_all50_07_snpEff_SIFT_formatted_${IDX}.vcf.gz* ${WORKDIR}/

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}
date "+%Y-%m-%d %T" >> ${LOG}

# ################################################################################

# ### Move intermediate files to Sirius
# cd ${WORKDIR}

# echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_all50_07_snpEff_SIFT_*.vcf.gz files to ${WORKSIRIDIR}... " >> ${PROGRESSLOG}

# scp JointCalls_all50_07_snpEff_SIFT_${IDX}.vcf.gz ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_all50_07_snpEff_${IDX}.vcf.gz FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# scp JointCalls_all50_07_snpEff_SIFT_${IDX}.vcf.gz.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_all50_07_snpEff_${IDX}.vcf.gz.tbi FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# # rm JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz
# # rm JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz.tbi
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done WGSproc7 for all50 ${REF}"
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done WGSproc7 for all50 ${REF}" >> ${PROGRESSLOG}


conda deactivate
