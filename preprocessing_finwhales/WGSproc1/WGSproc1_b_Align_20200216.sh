#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l h_rt=120:00:00,h_data=7G,h_vmem=15G,highp
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc1.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc1.err.txt
#$ -pe shared 6
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
NEXTSCRIPT=${SCRIPTDIR}/WGSproc1/WGSproc1_c_MergeAlign_Qualimap_20200216.sh  

# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc1_b_Align for ${NAME} ${REF} Job ID: ${JOB_ID}"

################################################################################

### get into working directory 

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}

PROGRESSLOG=WGSproc1_${RGID}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" >> ${PROGRESSLOG}

################################################################################

### AlignCleanBam
echo -e "[$(date "+%Y-%m-%d %T")] Aligning... " >> ${PROGRESSLOG}
LOG2=03_b_${RGID}_bwamem.log
date "+%Y-%m-%d %T" > ${LOG2}

## 2. bwa mem
bwa mem -M -t 6 -p -o ${RGID}_MIA_Aligned.bam \
${REFERENCE} ${RGID}_MarkIlluminaAdapters.fastq 2>> ${LOG2} 

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

rm ${RGID}_MarkIlluminaAdapters.fastq

date "+%Y-%m-%d %T" >> ${LOG2}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


################################################################################

### Submit next script if only one read set for sample (flag=0: do not qsub, flag=1: qsub)
# next script's input 
# 
# NAME=${1} # sample name
# RGID=${2} # Read group id 
# FLAG=${3} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
# USER=${4} # meixilin 
# REF=${5} # reference

if [ ${FLAG} -ne 0 ]; then
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc1_c_${NAME} ${NEXTSCRIPT} ${NAME} ${RGID} ${FLAG} ${USER} ${REF}) 
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Done WGSproc1_b_Align for ${NAME} ${REF} Job ID: ${JOB_ID}"

conda deactivate


