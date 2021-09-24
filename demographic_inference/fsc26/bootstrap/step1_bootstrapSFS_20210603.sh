#!/bin/bash
#$ -l highp,h_data=8G,h_vmem=10G,h_rt=23:59:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/fsc26/param_btsp/step1_bootstrapSFS_20210603.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/fsc26/param_btsp/step1_bootstrapSFS_20210603.err.txt
#$ -m abe

# @version 		v2
# @usage		qsub step1_param_bootstrap.sh <model> <submodel> <population> <prefix>
# @description	Generate bootstrapped SFS from given parameter estimations
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Feb 21 23:59:39 2021
# @modification: Tue Apr 27 10:13:48 2021
# @modification: Using new settings of nchr and nloci
# @modification: Thu Jun  3 20:05:14 2021
# @modification: Update based on the Total sites used
# Testing parameters:
# MODEL='1DModels'
# SUBMODEL='1D.1Epoch'
# POP='ENP'
# PREFIX='v1_r90'

############################################################
## import packages
sleep $((RANDOM % 120))

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
REFID=${1}
SETTING='neutral'

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
FSCFILE=${HOMEDIR}/scripts/fsc26/param_btsp/fsc_output_list.csv

if [[ ! -f ${FSCFILE} ]]; then
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: ${FSCFILE} does not exist"
    exit 1
fi

# get the parameter values directly from the reference ID
MODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $2}' $FSCFILE)
SUBMODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $3}' $FSCFILE)
POP=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $4}' $FSCFILE)
PREFIX=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $5}' $FSCFILE)
BOOTPREFIX=${SUBMODEL}.${POP}_${PREFIX}
NSIM=100
TODAY=$(date "+%Y%m%d")

WORKDIR=${HOMEDIR}/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}/${SUBMODEL}.${POP}/${BOOTPREFIX}_bootSFS
mkdir -p ${WORKDIR}
LOGDIR=${HOMEDIR}/reports/fsc26/param_btsp/step1_bootstrapSFS_20210603
mkdir -p ${LOGDIR}
LOG=${HOMEDIR}/reports/fsc26/param_btsp/step1_bootstrapSFS_20210603/step1_bootstrapSFS_${SETTING}_${BOOTPREFIX}_${TODAY}.log

# the bootstrap file
BOOTFILE=${HOMEDIR}/scripts/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}/${BOOTPREFIX}_boot.par
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# the software
FSC26=/u/project/rwayne/software/finwhale/fastsimcoal/fsc26_linux64/fsc26

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; Using ${BOOTFILE}; Outputting to ${WORKDIR}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; Using ${BOOTFILE}; Outputting to ${WORKDIR}" > ${LOG}

cd ${WORKDIR}

rsync -a ${BOOTFILE} ./

# -n  --numsims 1000      : number of simulations to perform. Also applies for parameter estimation
# -j  --jobs              : output one simulated or bootstrapped SFS per file
#                            in a separate directory for easier analysis
#                            (requires -d or -m and -s0 options)
# -m  --msfs              : computes minor site frequency spectrum
# -s  --dnatosnp 2000     : output DNA as SNP data, and specify maximum no. of SNPs
#                            to output (use 0 to output all SNPs).
# -x  --noarloutput       : does not generate Arlequin output
# -I  --inf               : generates DNA mutations according to an
#                            infinite site (IS) mutation model
# -q  --quiet             : minimal message output to console

${FSC26} -i ${BOOTPREFIX}_boot.par -n ${NSIM} -j -m -s0 -x -I -q &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}

