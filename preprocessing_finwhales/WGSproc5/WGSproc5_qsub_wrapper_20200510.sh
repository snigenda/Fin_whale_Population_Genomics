#!/bin/bash 

# Title: qsub wrapper for WGSproc5 Joint Genotyping 
# 
# Author: Meixi Lin
# Date: Thu Jan 30 12:20:44 2020

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc5/WGSproc5_GenotypeGVCFs_20200510.sh

USER=${1} ## "meixilin"/"snigenda"
REF=${2} # should be Minke all the time but for compatibility of the 6 individuals

if [ $REF == 'Minke' ]; then
    NJOBS=96
    HARD_RESOURCE="h_rt=23:00:00,h_data=23G"
fi
if [ $REF == 'Bryde' ]; then
    NJOBS=23
    HARD_RESOURCE="highp,h_rt=48:00:00,h_data=30G"
fi 

${QSUB} -t 1-${NJOBS} -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}

echo "The qsub script was"
echo "${QSUB} -t 1-${NJOBS} -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}"
