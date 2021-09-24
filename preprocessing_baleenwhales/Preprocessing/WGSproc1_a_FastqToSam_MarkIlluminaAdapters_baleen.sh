#! /bin/bash
#$ -l highp,h_rt=120:00:00,h_data=40G,h_vmem=45G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc1_20210112.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc1_20210112.err.txt
#$ -m abe

# @modification: Tue Jan 12 22:02:00 2021
# @modification: Update pipeline; stop backing up to sirius

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub


################################################################################

### Set variables

FQ1=${1} # forward read fastq.gz files, please use full path
FQ2=${2} # reverse read fastq.gz files, please use full path
NAME=${3} # for picard input: SAMPLE_NAME = Sample name to insert into the read group header Required.
RGID=${4} # for picard input: READ_GROUP_NAME = Read group name Default value: A.
RGLB=${5} # for picard input: LIBRARY_NAME = The library name to place into the LB attribute in the read group header
RGPU=${6} # for picard input: PLATFORM_UNIT = The platform unit (often run_barcode.lane) to insert into the read group header; {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_NAME}
RGCN=${7} # for picard input: SEQUENCING_CENTER = The sequencing center from which the data originated
RGPM=${8} # for picard input: PLATFORM_MODEL = "NovaSeq/HiSeq"
FLAG=${9} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
REF=${10} # reference name

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes
# SIRIUSDIR=/data3/finwhale/baleen_genomes

WORKDIR="${HOMEDIR}/preprocessing/${NAME}/WGSproc0_1a"

SCRIPTDIR="${HOMEDIR/baleen_genomes}/scripts/baleen_genomes"
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

NEXTSCRIPT=${SCRIPTDIR}/Preprocessing/WGSproc1_b_Align_baleen.sh

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Start WGSproc1a for ${NAME}; git commit id: ${COMMITID}"
echo "[$(date "+%Y-%m-%d %T")] The qsub input ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${REF}"

################################################################################

### Make directories

mkdir -p ${WORKDIR}
mkdir -p ${HOMEDIR}/preprocessing/${NAME}/${REF}
mkdir -p ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
# ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${SIRIUSDIR}/preprocessing"

cd ${SCRATCHDIR}/preprocessing/${NAME}
mkdir -p temp

PROGRESSLOG=${WORKDIR}/WGSproc1_a_FastqToSam_MarkIlluminaAdapters_baleen_${RGID}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}

################################################################################

### FastqToSam
echo -e "[$(date "+%Y-%m-%d %T")] Picard FastqToSam... " >> ${PROGRESSLOG}
LOG=${WORKDIR}/01_${RGID}_FastqToSam.log
date "+%Y-%m-%d %T" > ${LOG}

# CHECK: some variables not needed
picard -Xmx35G FastqToSam \
FASTQ=${FQ1} \
FASTQ2=${FQ2} \
OUTPUT=${RGID}_FastqToSam.bam \
READ_GROUP_NAME=${RGID} \
SAMPLE_NAME=${NAME} \
LIBRARY_NAME=${RGLB} \
PLATFORM_UNIT=${RGPU} \
SEQUENCING_CENTER=${RGCN} \
PLATFORM_MODEL=${RGPM} \
PLATFORM=illumina \
TMP_DIR=./temp 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


################################################################################

### MarkIlluminaAdapters
echo -e "[$(date "+%Y-%m-%d %T")] Picard MarkIlluminaAdapters... " >> ${PROGRESSLOG}
LOG=${WORKDIR}/02_${RGID}_MarkIlluminaAdapters.log
date "+%Y-%m-%d %T" > ${LOG}

picard -Xmx35G MarkIlluminaAdapters \
INPUT=${RGID}_FastqToSam.bam \
OUTPUT=${RGID}_MarkIlluminaAdapters.bam \
METRICS=${RGID}_MarkIlluminaAdapters_metrics.txt \
TMP_DIR=./temp 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

cp ${RGID}_MarkIlluminaAdapters_metrics.txt ${WORKDIR}/


################################################################################

### AlignCleanBam: Prepare

echo -e "[$(date "+%Y-%m-%d %T")] Prepare Aligning... " >> ${PROGRESSLOG}
LOG1=${WORKDIR}/03_a_${RGID}_SamToFastq.log
date "+%Y-%m-%d %T" > ${LOG1}

## 1. samtofastq
picard -Xmx35G SamToFastq \
INPUT=${RGID}_MarkIlluminaAdapters.bam \
FASTQ=${RGID}_MarkIlluminaAdapters.fastq \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=./temp 2>> ${LOG1}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
date "+%Y-%m-%d %T" >> ${LOG1}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

if [ ${FLAG} -ne 0 ]; then
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc1_b_${NAME}_${REF} ${NEXTSCRIPT} ${NAME} ${RGID} ${FLAG} ${REF})
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted WGSproc1_b_${NAME}_${REF} $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Done WGSproc1_a_FastqToSam_MarkIlluminaAdapters for ${NAME}"

conda deactivate


