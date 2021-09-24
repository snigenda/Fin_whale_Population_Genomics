#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -l h_rt=20:00:00,h_data=2G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc4.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc4.err.txt
#$ -m abe

set -o pipefail 
###########################################################

### Set variables
NAME=${1}
USER=${2}
REF=${3}
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

HOMEDIR=/u/project/rwayne/snigenda/finwhale
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale 
SIRIUSDIR=/data3/finwhale

if [ $REF == 'Minke' ]; then
    NJOBS=96
fi
if [ $REF == 'Bryde' ]; then
    NJOBS=23
fi 

cd ${SCRATCHDIR}/preprocessing/${NAME}/${REF}

PROGRESSLOG=WGSproc4_${NAME}_progress.log
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
echo -e "[$(date "+%Y-%m-%d %T")] Checking if all vcf have been moved to HOMEDIR... " >> ${PROGRESSLOG}

# Assign each number in the sequence to i; loop until we have done them all
for ((ii=1 ; ii<=${NJOBS} ; ii++)) ; do
   # Format the number with padding for the file name part
   printf -v id '%02d' "$ii"
   # If a file with this name does not exist, print error and exit 
   if [ ! -f "${HOMEDIR}/preprocessing/${NAME}/${REF}/${NAME}_${BAMHEAD}_${id}.g.vcf.gz" ] ; then
       # Print it to standard output
       echo -e "${HOMEDIR}/preprocessing/${NAME}/${REF}/${NAME}_${BAMHEAD}_${id}.g.vcf.gz is missing" >> ${PROGRESSLOG}
       echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
       exit 1
   fi
done
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

### Copy all intermediate logs to sirius 
echo -e "[$(date "+%Y-%m-%d %T")] Copying logs to Sirius... " >> ${PROGRESSLOG}
ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${SIRIUSDIR}/preprocessing/${NAME}/${REF}/scratchlogs"
scp *.log ${USER}@sirius.eeb.ucla.edu:${SIRIUSDIR}/preprocessing/${NAME}/${REF}/scratchlogs 

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

mkdir ${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs 
mv *log ${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs 
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi 

echo "${NAME},${REF},$(date "+%Y-%m-%d")" >> ${HOMEDIR}/preprocessing/preprocessing_donelist.txt

###########################################################

### Generate md5sum for the final folder 
cd ${HOMEDIR}/preprocessing/${NAME}/${REF}

PROGRESSLOG=${HOMEDIR}/preprocessing/${NAME}/${REF}/scratchlogs/WGSproc4_${NAME}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] Generating md5sum files ... " >> ${PROGRESSLOG}

find -type f -exec md5sum "{}" + > preprocessing_${NAME}_${REF}.md5sum
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi 

# remove the line that calculate the md5sum of md5sum 
sed -i '/.md5sum$/d' preprocessing_${NAME}_${REF}.md5sum
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi 

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

