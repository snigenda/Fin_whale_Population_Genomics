#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l h_rt=22:00:00,h_data=20G,h_vmem=24G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc1.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc1.err.txt
#$ -m abe

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
SIRIUSDIR=/data3/finwhale  

# we will decide which reference file to use for further analysis after testing mapping quality
if [ $REF == 'Minke' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
fi
if [ $REF == 'Bryde' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
fi

SCRIPTDIR=${HOMEDIR}/scripts
NEXTSCRIPT=${SCRIPTDIR}/WGSproc2/WGSproc2_MarkDuplicates_20200103.sh  

# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc1_c_MergeAlign_Qualimap for ${NAME} ${REF} Job ID: ${JOB_ID}"

################################################################################

### get into working directory 

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}

PROGRESSLOG=WGSproc1_${RGID}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" >> ${PROGRESSLOG}

################################################################################

### Merge Clean Bam
echo -e "[$(date "+%Y-%m-%d %T")] Merge Clean Bam... " >> ${PROGRESSLOG}
LOG3=03_c_${RGID}_MergeBamAlignment.log
date "+%Y-%m-%d %T" > ${LOG3}

## 3. mergebamalignment
picard -Xmx16G MergeBamAlignment \
ALIGNED_BAM=${RGID}_MIA_Aligned.bam \
UNMAPPED_BAM=${RGID}_FastqToSam.bam \
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
rm ${RGID}_FastqToSam.bam
rm ${RGID}_MarkIlluminaAdapters.bam

date "+%Y-%m-%d %T" >> ${LOG3}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

################################################################################

### Run qualimap on aligned bam
echo -e "[$(date "+%Y-%m-%d %T")] Qualimap... " >> ${PROGRESSLOG}
LOG=04_${RGID}_qualimap.log
date "+%Y-%m-%d %T" > ${LOG}

qualimap bamqc -bam ${RGID}_Aligned.bam -c --java-mem-size=16G

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

mv ${RGID}_Aligned_stats ${HOMEDIR}/preprocessing/${NAME}/${REF}
date "+%Y-%m-%d %T" >> ${LOG}
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
    NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc2_${NAME} ${NEXTSCRIPT} ${NAME} ${RGID} ${FLAG} ${USER} ${REF}) 
    echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}
fi

echo "[$(date "+%Y-%m-%d %T")] Done WGSproc1 for ${NAME} ${REF} Job ID: ${JOB_ID}"

conda deactivate


