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

NEXTSCRIPT=${SCRIPTDIR}/Preprocessing/WGSproc2_MarkDuplicates_baleen.sh

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Start WGSproc1_c_MergeAlign_Qualimap for ${NAME} ${REF}; git commit id: ${COMMITID}; The qsub input ${NAME} ${RGID} ${FLAG} ${REF}"

################################################################################

### Get into working directory

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}

PROGRESSLOG=WGSproc1_bc_Align_MergeQualimap_baleen_${RGID}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" >> ${PROGRESSLOG}

################################################################################

### Merge Clean Bam
echo -e "[$(date "+%Y-%m-%d %T")] Merge Clean Bam... " >> ${PROGRESSLOG}
LOG3=03_c_${RGID}_${REF}_MergeBamAlignment.log
date "+%Y-%m-%d %T" > ${LOG3}

## 3. mergebamalignment
picard -Xmx35G MergeBamAlignment \
ALIGNED_BAM=${RGID}_MIA_Aligned.bam \
UNMAPPED_BAM=${SCRATCHDIR}/preprocessing/${NAME}/${RGID}_FastqToSam.bam \
OUTPUT=${RGID}_Aligned.bam \
R=${REFERENCE} CREATE_INDEX=true \
ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=./temp 2>> ${LOG3}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

rm ${RGID}_MIA_Aligned.bam

date "+%Y-%m-%d %T" >> ${LOG3}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

################################################################################

### Run qualimap on aligned bam
echo -e "[$(date "+%Y-%m-%d %T")] Qualimap... " >> ${PROGRESSLOG}
LOG4=04_${RGID}_${REF}_qualimap.log
date "+%Y-%m-%d %T" > ${LOG4}

qualimap bamqc -bam ${RGID}_Aligned.bam -c --java-mem-size=35G &>> ${LOG4}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

mv ${RGID}_Aligned_stats ${HOMEDIR}/preprocessing/${NAME}/${REF}/
date "+%Y-%m-%d %T" >> ${LOG4}
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
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc2_${NAME}_${REF} ${NEXTSCRIPT} ${NAME} ${RGID} ${FLAG} ${REF})
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted WGSproc2_${NAME}_${REF} $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Done WGSproc1_c_MergeAlign_Qualimap for ${NAME} ${REF}"

conda deactivate


