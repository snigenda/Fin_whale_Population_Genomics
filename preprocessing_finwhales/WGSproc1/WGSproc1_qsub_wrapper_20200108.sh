#!/bin/bash 

# Title: qsub wrapper for process 1
# 
# Author: Meixi Lin
# Date: Tue Nov 19 14:53:25 2019

# FQ1=${1} # forward read fastq.gz files, please use full path
# FQ2=${2} # reverse read fastq.gz files, please use full path
# NAME=${3} # for picard input: SAMPLE_NAME = Sample name to insert into the read group header Required.
# RGID=${4} # for picard input: READ_GROUP_NAME = Read group name Default value: A.
# RGLB=${5} # for picard input: LIBRARY_NAME = The library name to place into the LB attribute in the read group header
# RGPU=${6} # for picard input: PLATFORM_UNIT = The platform unit (often run_barcode.lane) to insert into the read group header; {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_NAME}
# RGCN=${7} # for picard input: SEQUENCING_CENTER = The sequencing center from which the data originated 
# RGPM=${8} # for picard input: PLATFORM_MODEL = "NovaSeq/HiSeq"
# FLAG=${9} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
# USER=${10} # hoffman2 user name 
# REF=${11} # reference genome name, takes 'Minke'/'Bryde'

# 1. set variables for this environment ########

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/snigenda/finwhale
DICT=${HOMEDIR}/scripts/config/fqpath_fqname_rgid.csv # location of file
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc1/WGSproc1_a_FastqToSam_MarkIlluminaAdapters_20200216.sh

NAME=${1} ## only need to change this line for submitting different individual's file
FLAG=${2} ## 0 if only run this section, 1 if keep going 
USER=${3} ## "meixilin"/"snigenda"
REF=${4} ## "Minke"/

# 2. Automatic, no change needed: generate input for the worker script ########
# explanation on awk: the $DICT is a csv that contains quoted columns
# it searches in the second column for match and prints the corresponding other columns value
# "AccessionID" "SampleName"  "R1"          "R2"          "RGID"        "RGLB"        "RGPU"        "RGCN"        "RGPM"
FQ1=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $3}' $DICT | tr -d \"`
FQ2=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $4}' $DICT | tr -d \"`
RGID=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $5}' $DICT | tr -d \"` 
RGLB=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $6}' $DICT | tr -d \"` 
RGPU=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $7}' $DICT | tr -d \"` 
RGCN=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $8}' $DICT | tr -d \"`
RGPM=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $9}' $DICT | tr -d \"`


$QSUB -M $USER $WORKSCRIPT ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${USER} ${REF}

# echo the input 
echo "The qsub script was"
echo "$QSUB -M $USER $WORKSCRIPT ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${USER} ${REF}"


