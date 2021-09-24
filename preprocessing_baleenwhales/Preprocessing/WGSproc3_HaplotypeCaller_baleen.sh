#! /bin/bash
#$ -l highp,h_rt=72:00:00,h_data=20G,h_vmem=24G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc3_20210112.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc3_20210112.err.txt
#$ -m abe

# @version         v1
# @usage           qsub -t 1-${NJOBS} -N WGSproc3_${NAME}_${REF} WGSproc3_HaplotypeCaller_baleen.sh ${NAME} ${REF}
# @description     Generate the haplotypes
# @modification: Tue Jan 12 23:54:26 2021
# @modification: Update pipeline; stop backing up to sirius

sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

# QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub

################################################################################

### Set variables
NAME=${1}
REF=${2}
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes
# SIRIUSDIR=/data3/finwhale/baleen_genomes

IDX=$(printf %02d ${SGE_TASK_ID})

if [ $REF == 'Minke' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${REFERENCE/GCF_000493695.1_BalAcu1.0_genomic.fasta/contiglist\/BalAcu1.0_genomic.contiglist}_${IDX}.list
fi
if [ $REF == 'Bryde' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
    INTERVAL=${REFERENCE/Balaenoptera_edeni_HiC.fasta/contiglist\/Bal_edeni_HiC.contiglist}_${IDX}.list
fi
if [ $REF == 'Humpback' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/humpback_whale_genome/GCA_004329385.1_megNov1/GCA_004329385.1_megNov1_genomic.fasta
    INTERVAL=${REFERENCE/GCA_004329385.1_megNov1_genomic.fasta/contiglist\/GCA_004329385.1_megNov1_genomic.contiglist}_${IDX}.list
fi
if [ $REF == 'Blue' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta
    INTERVAL=${REFERENCE/GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta/contiglist\/GCF_009873245.2_mBalMus1.pri.v3_genomic.contiglist}_${IDX}.list
fi

SCRIPTDIR="${HOMEDIR/baleen_genomes}/scripts/baleen_genomes"
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

################################################################################

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc3_HaplotypeCaller for ${NAME} ${REF}; git commit id: ${COMMITID}; The qsub input ${NAME} ${REF}"

mkdir -p ${HOMEDIR}/preprocessing/${NAME}/${REF}
mkdir -p ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
# ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${SIRIUSDIR}/preprocessing/${NAME}/${REF}"
cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
mkdir -p temp

PROGRESSLOG=WGSproc3_HaplotypeCaller_${NAME}_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${IDX}" > ${PROGRESSLOG}

################################################################################

### Check to make sure that alignment jobs completed successfully
echo -e "[$(date "+%Y-%m-%d %T")] Checking previous job completion status... " >> ${PROGRESSLOG}

FCOUNT=$(cat ${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs/*progress.log | grep "FAIL" | wc -l)

if [ $FCOUNT -ne 0 ]
then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi


################################################################################

### Generate gVCF files
echo -e "[$(date "+%Y-%m-%d %T")] GATK HaplotypeCaller... " >> ${PROGRESSLOG}
LOG=06_${NAME}_${REF}_${BAMHEAD}_HaplotypeCaller_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx15G -Djava.io.tmpdir=./temp -T HaplotypeCaller \
-R ${REFERENCE} \
-ERC BP_RESOLUTION \
-mbq 20 \
-out_mode EMIT_ALL_SITES \
-I ${HOMEDIR}/preprocessing/${NAME}/${REF}/${NAME}_${BAMHEAD}.bam \
-L ${INTERVAL} \
-o ${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


### Save gVCF files
echo -e "[$(date "+%Y-%m-%d %T")] Moving gVCF files to ${HOMEDIR}... " >> ${PROGRESSLOG}

# scp ${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz* ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}/
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] scp FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

mv ${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz* ${HOMEDIR}/preprocessing/${NAME}/${REF}/
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] mv FAIL" >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc3_HaplotypeCaller for ${NAME} ${REF}"

conda deactivate


