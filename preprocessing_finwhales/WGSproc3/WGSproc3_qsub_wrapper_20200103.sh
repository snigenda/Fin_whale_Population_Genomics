#!/bin/bash 

# Title: qsub wrapper for WGSproc3 HaplotypeCaller on ENPAK28
# 
# Author: Meixi Lin
# Date: Tue Dec 10 16:10:29 PST 2019

# ### Set variables
# NAME={1}
# USER={2}
# REF="Minke"
# BAMHEAD="RemoveBadReads" # MarkDuplicate/RemoveBadReads

# 1. set variables for this environment ########

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc3/WGSproc3_HaplotypeCaller_20200103.sh

NAME=${1} ## only need to change this line for submitting different individual's file
USER=${2} ## "meixilin"/"snigenda"
REF=${3} # should be Minke all the time but for compatibility of the 6 individuals

if [ $REF == 'Minke' ]; then
    NJOBS=96
fi
if [ $REF == 'Bryde' ]; then
    NJOBS=23
fi 

${QSUB} -t 1-${NJOBS} ${WORKSCRIPT} ${NAME} ${USER} ${REF}

echo "The qsub script was"
echo "${QSUB} -t 1-${NJOBS} ${WORKSCRIPT} ${NAME} ${USER} ${REF}"