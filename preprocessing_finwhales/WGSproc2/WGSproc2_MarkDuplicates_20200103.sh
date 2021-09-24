#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l h_rt=22:00:00,h_data=10G,h_vmem=16G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc2.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc2.err.txt
#$ -m abe

# this step: marks duplicate and remove bad reads 
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub


################################################################################

### Set variables

NAME=${1} # sample name
RGID=${2} # Read group id 
FLAG=${3} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
USER=${4} # meixilin 
REF=${5} # reference


HOMEDIR=/u/project/rwayne/snigenda/finwhale
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale 
SIRIUSDIR=/data3/finwhale  ## CHECK: we will need to check the location and name of this directory and of the files inside

# we will decide which reference file to use for further analysis after testing mapping quality
if [ $REF == 'Minke' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    NJOBS=96
fi
if [ $REF == 'Bryde' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
    NJOBS=23
fi

SCRIPTDIR=${HOMEDIR}/scripts
NEXTSCRIPT=${SCRIPTDIR}/WGSproc3/WGSproc3_HaplotypeCaller_20200103.sh 

# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc2 for ${NAME} ${REF} Job ID: ${JOB_ID}"
echo "The qsub input"
echo "${NAME} ${RGID} ${FLAG} ${USER} ${REF}"

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF} 
mkdir -p temp

PROGRESSLOG=WGSproc2_${NAME}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


################################################################################

### Check to make sure that alignment jobs completed successfully
echo -e "[$(date "+%Y-%m-%d %T")] Checking alignment job completion status... " >> ${PROGRESSLOG}

FCOUNT=$(cat *progress.log | grep "FAIL" | wc -l)

if [ $FCOUNT -ne 0 ]
then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi


################################################################################

### MarkDuplicates
echo -e "[$(date "+%Y-%m-%d %T")] Picard MarkDuplicates... " >> ${PROGRESSLOG}
LOG=05_a_${NAME}_MarkDuplicates.log
date "+%Y-%m-%d %T" > ${LOG}

picard -Xmx8G -Djava.io.tmpdir=./temp MarkDuplicates \
INPUT=${RGID}_Aligned.bam \
OUTPUT=${NAME}_MarkDuplicates.bam \
METRICS_FILE=${NAME}_MarkDuplicates_metrics.txt \
MAX_RECORDS_IN_RAM=150000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
CREATE_INDEX=true \
TMP_DIR=./temp 2>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mv ${NAME}_MarkDuplicates_metrics.txt ${HOMEDIR}/preprocessing/${NAME}/${REF}


### Copying MarkDuplicates.bam to Sirius and move to HOMEDIR
echo -e "[$(date "+%Y-%m-%d %T")] Copying MarkDuplicates.bam to Sirius/HOMEDIR... " >> ${PROGRESSLOG}

scp ${NAME}_MarkDuplicates.bam ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}
exit1=${?}
scp ${NAME}_MarkDuplicates.bai ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}
exit2=${?}
mv ${NAME}_MarkDuplicates.bam ${HOMEDIR}/preprocessing/${NAME}/${REF}
exit3=${?}
mv ${NAME}_MarkDuplicates.bai ${HOMEDIR}/preprocessing/${NAME}/${REF}
exit4=${?}

let "exitVal=$exit1+$exit2+$exit3+$exit4" 
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


rm ${RGID}_Aligned.bam
rm ${RGID}_Aligned.bai


################################################################################

### Submit next scripts

if [ ${FLAG} -ne 0 ]; then
    NEXT_JOB_ID=$(${QSUB} -terse -t 1-${NJOBS} -N WGSproc3_${NAME} ${NEXTSCRIPT} ${NAME} ${USER} ${REF} | cut -d'.' -f1)
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Done WGSproc2 for ${NAME} ${REF} Job ID: ${JOB_ID}"

conda deactivate
