#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l highp,h_rt=48:00:00,h_data=23G,h_vmem=23G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc3.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc3.err.txt
#$ -m abe


source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

# QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub

################################################################################

### Set variables
NAME=${1}
USER=${2}
REF=${3}
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

HOMEDIR=/u/project/rwayne/snigenda/finwhale
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale 
SIRIUSDIR=/data3/finwhale 
mkdir -p ${SCRATCHDIR}/preprocessing/${NAME}/${REF}

IDX=$(printf %02d ${SGE_TASK_ID})


if [ $REF == 'Minke' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${REFERENCE/GCF_000493695.1_BalAcu1.0_genomic.fasta/contiglist\/BalAcu1.0_genomic.contiglist}_${IDX}.list
fi

if [ $REF == 'Bryde' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
    INTERVAL=${REFERENCE/Balaenoptera_edeni_HiC.fasta/contiglist\/Bal_edeni_HiC.contiglist}_${IDX}.list
fi

# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc3 for ${NAME} ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"
echo "The qsub input"
echo "${NAME} ${USER} ${REF} ${SGE_TASK_ID}"


cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}
mkdir -p temp

PROGRESSLOG=WGSproc3_${NAME}_${REF}_${BAMHEAD}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}

################################################################################

### Check to make sure that alignment jobs completed successfully
echo -e "[$(date "+%Y-%m-%d %T")] Checking previous job completion status... " >> ${PROGRESSLOG}

FCOUNT=$(cat *progress.log | grep "FAIL" | wc -l)

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

gatk3 -Xmx15g -Djava.io.tmpdir=./temp -T HaplotypeCaller \
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
echo -e "[$(date "+%Y-%m-%d %T")] Moving gVCF files to sirius/${HOMEDIR}... " >> ${PROGRESSLOG}

scp ${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz* ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] scp FAIL" >> ${PROGRESSLOG}
    exit 1
fi
mv ${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz* ${HOMEDIR}/preprocessing/${NAME}/${REF}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] mv FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


echo "[$(date "+%Y-%m-%d %T")] Done WGSproc3 for ${NAME} ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"
conda deactivate 


