#! /bin/bash
#$ -l highp,h_rt=120:00:00,h_data=40G,h_vmem=45G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc2_20210112.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc2_20210112.err.txt
#$ -m abe

# @modification: Tue Jan 12 22:02:00 2021
# @modification: Update pipeline; stop backing up to sirius

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
REF=${4} # reference

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes
# SIRIUSDIR=/data3/finwhale/baleen_genomes

if [ $REF == 'Minke' ]; then
    NJOBS=96
    REFERENCE=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
fi
if [ $REF == 'Bryde' ]; then
    NJOBS=23
    REFERENCE=${REFDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
fi
if [ $REF == 'Humpback' ]; then
    NJOBS=73
    REFERENCE=${REFDIR}/cetacean_genomes/humpback_whale_genome/GCA_004329385.1_megNov1/GCA_004329385.1_megNov1_genomic.fasta
fi
if [ $REF == 'Blue' ]; then
    NJOBS=24
    REFERENCE=${REFDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta
fi

SCRIPTDIR="${HOMEDIR/baleen_genomes}/scripts/baleen_genomes"
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

NEXTSCRIPT=${SCRIPTDIR}/Preprocessing/WGSproc3_HaplotypeCaller_baleen.sh

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Start WGSproc2_MarkDuplicates for ${NAME} ${REF}; git commit id: ${COMMITID}; The qsub input ${NAME} ${RGID} ${FLAG} ${REF}"

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
mkdir -p temp

PROGRESSLOG=WGSproc2_MarkDuplicates_${NAME}_${REF}_progress.log
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
LOG=05_a_${NAME}_${REF}_MarkDuplicates.log
date "+%Y-%m-%d %T" > ${LOG}

picard -Xmx35G -Djava.io.tmpdir=./temp MarkDuplicates \
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

cp ${NAME}_MarkDuplicates_metrics.txt ${HOMEDIR}/preprocessing/${NAME}/${REF}/


### Moving MarkDuplicates.bam to HOMEDIR
echo -e "[$(date "+%Y-%m-%d %T")] Moving MarkDuplicates.bam to ${HOMEDIR}... " >> ${PROGRESSLOG}

# scp ${NAME}_MarkDuplicates.bam ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}
# exit1=${?}
# scp ${NAME}_MarkDuplicates.bai ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}
# exit2=${?}
mv ${NAME}_MarkDuplicates.bam ${HOMEDIR}/preprocessing/${NAME}/${REF}/
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"  >> ${PROGRESSLOG}
    exit 1
fi

mv ${NAME}_MarkDuplicates.bai ${HOMEDIR}/preprocessing/${NAME}/${REF}/
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"  >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

###########################################################

### Copying logs
mkdir -p ${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs
cp *.log ${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs/

rm ${RGID}_Aligned.bam
rm ${RGID}_Aligned.bai


################################################################################

### Submit next scripts

if [ ${FLAG} -ne 0 ]; then
    NEXT_JOB_ID=$(${QSUB} -terse -t 1-${NJOBS} -N WGSproc3_${NAME}_${REF} ${NEXTSCRIPT} ${NAME} ${REF} | cut -d'.' -f1)
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted WGSproc3_${NAME}_${REF} $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Done WGSproc2_MarkDuplicates for ${NAME} ${REF}"

conda deactivate
