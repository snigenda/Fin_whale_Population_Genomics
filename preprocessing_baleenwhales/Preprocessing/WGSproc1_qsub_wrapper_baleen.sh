#!/bin/bash
#
# @version 		v1
# @script		WGSproc1_qsub_wrapper_baleen.sh
# @description	Wrapper for the baleen comparison preprocessing
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Aug  3 01:00:33 2020
# @modification: Tue Jan 12 22:09:23 2021
# @modification: Update pipeline; remove predicted size

###########################################################
## import packages

###########################################################
## def functions

###########################################################
## def variables
# FQ1=${1} # forward read fastq.gz files, please use full path
# FQ2=${2} # reverse read fastq.gz files, please use full path
# NAME=${3} # for picard input: SAMPLE_NAME = Sample name to insert into the read group header Required.
# RGID=${4} # for picard input: READ_GROUP_NAME = Read group name Default value: A.
# RGLB=${5} # for picard input: LIBRARY_NAME = The library name to place into the LB attribute in the read group header
# RGPU=${6} # for picard input: PLATFORM_UNIT = The platform unit (often run_barcode.lane) to insert into the read group header; {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_NAME}
# RGCN=${7} # for picard input: SEQUENCING_CENTER = The sequencing center from which the data originated
# RGPM=${8} # for picard input: PLATFORM_MODEL = "NovaSeq/HiSeq"
# FLAG=${9} # 0/1 whether or not to continue the processing (if 0, don't continue; if 1, continue)
# REF=${10} # reference name

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
DICT=${HOMEDIR}/scripts/config/baleen_fqpath_fqname.csv # location of file
WORKSCRIPT=${HOMEDIR}/scripts/baleen_genomes/Preprocessing/WGSproc1_a_FastqToSam_MarkIlluminaAdapters_baleen.sh

NAME=${1} ## only need to change this line for submitting different individual's file
FLAG=${2} ## 0 if only run this section, 1 if keep going
REF=${3} ## Minke/Bryde/Humpback/Blue

# 2. Automatic, no change needed: generate input for the worker script ########
# explanation on awk: the $DICT is a csv that contains quoted columns
# it searches in the second column for match and prints the corresponding other columns value
# SRA_Accession,SampleName,R1,R2,RGID,RGLB,RGPU,RGCN,RGPM,Organism,NCBISampleName,AvgSpotLen,Bases,ftpR1,ftpR2

FQ1=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $3}' $DICT)
FQ2=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $4}' $DICT)
RGID=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $5}' $DICT)
RGLB=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $6}' $DICT)
RGPU=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $7}' $DICT)
RGCN=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $8}' $DICT)
RGPM=$(awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $9}' $DICT)

###########################################################
## main
$QSUB -N WGSproc1_a_${NAME}_${REF} $WORKSCRIPT ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${REF}

# echo the input
echo "The qsub script was"
echo "$QSUB $WORKSCRIPT ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${REF}"
