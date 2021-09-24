#! /bin/bash
#$ -l highp,h_rt=120:00:00,h_data=7G,h_vmem=8G
#$ -pe shared 6
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

NAME=${1} # sample name
RGID=${2} # Read group id
FLAG=${3} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
REF=${4} # reference

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes
# SIRIUSDIR=/data3/finwhale/baleen_genomes

if [ $REF == 'Minke' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
fi
if [ $REF == 'Bryde' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
fi
if [ $REF == 'Humpback' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/humpback_whale_genome/GCA_004329385.1_megNov1/GCA_004329385.1_megNov1_genomic.fasta
fi
if [ $REF == 'Blue' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta
fi


SCRIPTDIR="${HOMEDIR/baleen_genomes}/scripts/baleen_genomes"
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

NEXTSCRIPT=${SCRIPTDIR}/Preprocessing/WGSproc1_c_MergeAlign_Qualimap_baleen.sh

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Start WGSproc1_b_Align for ${NAME} ${REF}; git commit id: ${COMMITID}; The qsub input ${NAME} ${RGID} ${FLAG} ${REF}"

################################################################################

### Get into working directory

mkdir -p ${HOMEDIR}/preprocessing/${NAME}/${REF}
mkdir -p ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
# ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${SIRIUSDIR}/preprocessing/${NAME}/${REF}"

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
mkdir -p temp

PROGRESSLOG=WGSproc1_bc_Align_MergeQualimap_baleen_${RGID}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}

################################################################################

### AlignCleanBam
echo -e "[$(date "+%Y-%m-%d %T")] Aligning... " >> ${PROGRESSLOG}
LOG2=03_b_${RGID}_${REF}_bwamem.log
date "+%Y-%m-%d %T" > ${LOG2}

## 2. bwa mem
bwa mem -M -t 6 -p -o ${RGID}_MIA_Aligned.bam \
${REFERENCE} ../${RGID}_MarkIlluminaAdapters.fastq 2>> ${LOG2}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

# rm ${RGID}_MarkIlluminaAdapters.fastq

date "+%Y-%m-%d %T" >> ${LOG2}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


################################################################################

### Submit next script if only one read set for sample (flag=0: do not qsub, flag=1: qsub)
# next script's input
#
# NAME=${1} # sample name
# RGID=${2} # Read group id
# FLAG=${3} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
# REF=${4} # reference

if [ ${FLAG} -ne 0 ]; then
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc1_c_${NAME}_${REF} ${NEXTSCRIPT} ${NAME} ${RGID} ${FLAG} ${REF})
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted WGSproc1_c_${NAME}_${REF} $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Done WGSproc1_b_Align for ${NAME} ${REF}"

conda deactivate


