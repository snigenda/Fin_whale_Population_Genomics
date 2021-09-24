#! /bin/bash
#$ -l h_rt=20:00:00,h_data=4G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc4_20210119.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc4_20210119.err.txt
#$ -m abe

# @version      v0
# @script       qsub WGSproc4_FinalCheck_baleen.sh <name> <ref>
# @description  Clean up the files and logs
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Jan 18 13:26:50 2021

###########################################################
## import packages
set -o pipefail

###########################################################
## def functions

###########################################################
## def variables
NAME=${1}
REF=${2}
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes
LOGDIR=${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs

if [ $REF == 'Minke' ]; then
    NJOBS="96"
    REFERENCE=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
fi
if [ $REF == 'Bryde' ]; then
    NJOBS="23"
    REFERENCE=${REFDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
fi

SCRIPTDIR="${HOMEDIR/baleen_genomes}/scripts/baleen_genomes"
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

###########################################################
## main

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}

PROGRESSLOG=${LOGDIR}/WGSproc4_FinalCheck_${NAME}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Start WGSproc4_FinalCheck for ${NAME} ${REF}; git commit id: ${COMMITID}; The qsub input ${NAME} ${REF}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}

################################################################################

### Check to make sure that all jobs completed successfully
echo -e "[$(date "+%Y-%m-%d %T")] Checking completion status of all jobs... " >> ${PROGRESSLOG}

FCOUNT=$(cat *progress.log | grep "FAIL" | wc -l)

if [ $FCOUNT -ne 0 ]
then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

################################################################################

### Check to make sure that all vcf have been moved to $HOMEDIR

echo -e "[$(date "+%Y-%m-%d %T")] Checking if all vcf have been moved to ${HOMEDIR}... " >> ${PROGRESSLOG}

# Loop through each contigs
for ((ii=1 ; ii<=${NJOBS} ; ii++)) ; do
    # Format the number with padding for the file name part
    IDX=$(printf %02d ${ii})
    # If a file with this name does not exist, print error and exit
    if [ ! -f "${HOMEDIR}/preprocessing/${NAME}/${REF}/${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz" ]
    then
        # Print it to standard output
        echo -e "${HOMEDIR}/preprocessing/${NAME}/${REF}/${NAME}_${BAMHEAD}_${IDX}.g.vcf.gz is missing" >> ${PROGRESSLOG}
        echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
        exit 1
    fi
done

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

### Copy all intermediate logs to $HOMEDIR
echo -e "[$(date "+%Y-%m-%d %T")] Archiving log files to ${LOGDIR}... " >> ${PROGRESSLOG}
rsync -ah --checksum --update *log ${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs/
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo "${NAME},${REF},$(date "+%Y-%m-%d")" >> ${HOMEDIR}/preprocessing/preprocessing_donelist.txt

###########################################################

### Generate md5sum for the final folder
cd ${HOMEDIR}/preprocessing/${NAME}/${REF}

echo -e "[$(date "+%Y-%m-%d %T")] Generating md5sum files... " >> ${PROGRESSLOG}

find -type f -exec md5sum "{}" + > preprocessing_${NAME}_${REF}.md5sum
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

# remove the line that calculate the md5sum of md5sum
sed -i '/.md5sum$/d' preprocessing_${NAME}_${REF}.md5sum
sed -i '/WGSproc4_FinalCheck/d' preprocessing_${NAME}_${REF}.md5sum

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Done WGSproc4_FinalCheck for ${NAME} ${REF}"

