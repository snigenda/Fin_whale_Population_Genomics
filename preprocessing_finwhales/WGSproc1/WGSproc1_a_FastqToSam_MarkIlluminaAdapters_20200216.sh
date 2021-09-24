#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l h_rt=22:00:00,h_data=24G,h_vmem=30G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc1.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc1.err.txt
#$ -m abe

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
USER=${10} # hoffman2 user name 
REF=${11} # reference genome name, takes 'Minke'/'Bryde'

HOMEDIR=/u/project/rwayne/snigenda/finwhale
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale 
SIRIUSDIR=/data3/finwhale  

# we will decide which reference file to use for further analysis after testing mapping quality
if [ $REF == 'Minke' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
fi
if [ $REF == 'Bryde' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
fi

SCRIPTDIR=${HOMEDIR}/scripts
NEXTSCRIPT=${SCRIPTDIR}/WGSproc1/WGSproc1_b_Align_20200216.sh  

# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc1 for ${NAME} ${REF} Job ID: ${JOB_ID}"
echo "The qsub input"
echo "${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${USER} ${REF}"

################################################################################

### Make directories

mkdir -p ${HOMEDIR}/preprocessing
mkdir -p ${SCRATCHDIR}/preprocessing

mkdir -p ${HOMEDIR}/preprocessing/${NAME}/${REF}
mkdir -p ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
# chmod -R 775 ${HOMEDIR}/preprocessing/${NAME}/${REF} # set permissions

ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${SIRIUSDIR}/preprocessing" 
ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${SIRIUSDIR}/preprocessing/${NAME}/${REF}" 
# ssh ${USER}@sirius.eeb.ucla.edu "chmod -R 777 ${SIRIUSDIR}/preprocessing/${NAME}/${REF}" 


cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
mkdir -p temp

PROGRESSLOG=WGSproc1_${RGID}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}

################################################################################

### FastqToSam
echo -e "[$(date "+%Y-%m-%d %T")] Picard FastqToSam... " >> ${PROGRESSLOG}
LOG=01_${RGID}_FastqToSam.log
date "+%Y-%m-%d %T" > ${LOG}

# -Xmx40G limits java's memory usage
# CHECK: some variables not needed 
picard -Xmx16G FastqToSam \
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
LOG=02_${RGID}_MarkIlluminaAdapters.log
date "+%Y-%m-%d %T" > ${LOG}

picard -Xmx16G MarkIlluminaAdapters \
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

mv ${RGID}_MarkIlluminaAdapters_metrics.txt ${HOMEDIR}/preprocessing/${NAME}/${REF}


################################################################################

### AlignCleanBam: Prepare

echo -e "[$(date "+%Y-%m-%d %T")] Prepare Aligning... " >> ${PROGRESSLOG}
LOG1=03_a_${RGID}_SamToFastq.log
date "+%Y-%m-%d %T" > ${LOG1}

## 1. samtofastq
picard -Xmx16G SamToFastq \
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
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc1_b_${NAME} ${NEXTSCRIPT} ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${USER} ${REF}) ## we wil check if this nomenclature it is maintained in our case
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Done WGSproc1_a_Prepare_Align for ${NAME} ${REF} Job ID: ${JOB_ID}"

conda deactivate


