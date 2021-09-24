#!/bin/bash 

# Title: qsub wrapper for process 2
# 
# Author: Meixi Lin
# Date: Mon Dec  2 00:03:17 PST 2019

# NAME=${1} # sample name
# RGID=${2} # Read group id 
# FLAG=${3} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
# USER=${4} # meixilin 
# REF=${5} # reference

# 1. set variables for this environment ########

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/snigenda/finwhale
DICT=${HOMEDIR}/scripts/config/fqpath_fqname_rgid.csv # location of file
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc2/WGSproc2_MarkDuplicates_20200103.sh

NAME=${1} ## only need to change this line for submitting different individual's file
FLAG=${2} ## 0 if only run this section, 1 if keep going 
USER=${3} ## "meixilin"/"snigenda"
REF=${4} ## "Minke"/ "Bryde"

# 2. Automatic, no change needed: generate input for the worker script ########
# explanation on awk: the $DICT is a csv that contains quoted columns
# it searches in the second column for match and prints the corresponding other columns value
# "AccessionID" "SampleName"  "R1"          "R2"          "RGID"        "RGLB"        "RGPU"        "RGCN"        "RGPM"
RGID=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $5}' $DICT | tr -d \"` 

$QSUB -M $USER $WORKSCRIPT ${NAME} ${RGID} ${FLAG} ${USER} ${REF}

# echo the input 
echo "The qsub script was"
echo "$QSUB -M $USER $WORKSCRIPT ${NAME} ${RGID} ${FLAG} ${USER} ${REF}"


