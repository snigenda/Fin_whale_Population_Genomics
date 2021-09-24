#!/bin/bash
#$ -l h_data=8G,h_rt=23:00:00
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc9_b_get_md5sum_all50_20201206.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc9_b_get_md5sum_all50_20201206.err.txt
#$ -m abe

# @version 		v2
# @usage		get the md5sum summary of the file
# @description	WGSproc9_all50
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Aug 14 10:58:02 2020
# @modification Mon Aug 24 23:04:20 2020
# @modification Final archive
# @modification: Sun Dec  6 17:09:06 2020
# @modification: Update to match all50 November 2020 version

###########################################################
## import packages
set -o pipefail

###########################################################
## def functions

###########################################################
## def variables
REF="Minke"
TODAY=$(date "+%Y%m%d")
LOCAL=/u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke/
OUTPREFIX="filteredvcf_all50_Minke_${TODAY}"

###########################################################
## main
cd ${LOCAL}

# check previous job completion
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Checking previous job completion status... "
FCOUNT=$(cat logs/*progress_all50.log | grep "FAIL" | wc -l)
if [ $FCOUNT -ne 0 ]
then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] SUCCESS previous job completion status... "

# get md5sum
echo -e "[$(date "+%Y-%m-%d %T")] Gathering md5sum file... "
# get md5sum for hoffman2 files (not including the md5sum file itself)
find -type f -exec md5sum "{}" + > ${OUTPREFIX}.md5sum
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
sed -i "/${OUTPREFIX}.md5sum$/d" ${OUTPREFIX}.md5sum
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] SUCCESS Gathering md5sum file... " >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

