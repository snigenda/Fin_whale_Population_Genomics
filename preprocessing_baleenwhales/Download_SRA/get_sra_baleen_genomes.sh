#!/bin/bash
#$ -l highp,h_data=8G,h_rt=120:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/get_sra_baleen_genomes_20210108.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/get_sra_baleen_genomes_20210108.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub get_sra_baleen_genomes.sh <SRAID>
# @description	get_sra_baleen_genomes.sh
# @modification Fri Jul 10 10:00:07 2020
# @modification change the srapath step to wget
# @modification: Fri Jan  8 16:45:52 2021
# @modification: 1. add md5sum check for SRA and fastq files 2. move to ENA system
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jun 26 12:27:38 2020

###########################################################
## import packages
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
SRAID=${1}

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/sra_seq
LOGDIR=/u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes

LOG=${LOGDIR}/get_sra_baleen_genomes_${SRAID}.log
MD5FILE=${LOGDIR}/get_sra_baleen_genomes.md5sum
DICT=${HOMEDIR}/scripts/config/baleen_fqpath_fqname.csv

SRAPATH_R1=$(awk -v pat="$SRAID" 'BEGIN {FS = ","}; $1 ~ pat {print $14}' $DICT)
SRAPATH_R2=$(awk -v pat="$SRAID" 'BEGIN {FS = ","}; $1 ~ pat {print $15}' $DICT)

COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
###########################################################
## main
mkdir -p $SCRATCHDIR
cd $SCRATCHDIR

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; SRA ID ${SRAID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; SRA ID ${SRAID}" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] wget ${SRAPATH_R1}" >> ${LOG}
wget -nv ${SRAPATH_R1} &>> ${LOG}
exitVal1=${?}
echo -e "[$(date "+%Y-%m-%d %T")] Done wget ${SRAPATH_R1}, exitVal = ${exitVal1}" >> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] wget ${SRAPATH_R2}" >> ${LOG}
wget -nv ${SRAPATH_R2} &>> ${LOG}
exitVal2=${?}
echo -e "[$(date "+%Y-%m-%d %T")] Done wget ${SRAPATH_R2}, exitVal = ${exitVal2}" >> ${LOG}

# get the md5sum
md5sum $SCRATCHDIR/${SRAID}_[0-9].fastq.gz >> ${MD5FILE}

echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}

conda deactivate