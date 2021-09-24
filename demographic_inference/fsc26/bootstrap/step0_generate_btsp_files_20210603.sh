#!/bin/bash
#
# @version      v0
# @script       bash step0_generate_btsp_files_20210428.sh <REMOTE/LOCAL> <SETTING> <DIR> <PREFIX> <VERSION> <RUN>
# @description  generate bootstrap files for fsc26 bootstraps
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Apr 28 13:02:09 2021
# Testing:
# LOCAL='local'
# SETTING='neutral'
# REFID='1'

###########################################################
## import packages
set -xeuo pipefail

###########################################################
## def functions

###########################################################
## def variables
LOCAL=${1} # remote/local
SETTING=${2} # neutral/whole
REFID=${3} # which line from fsc_output_list.csv to read in

if [ "${LOCAL}" == 'local' ]; then
    SOURCE0=/Users/linmeixi/google_drive/finwhale/analyses/fsc26/output
    SINK0=/Users/linmeixi/Lab/fin_whale/scripts_analyses/fsc26/param_btsp/${SETTING}
    FSCFILE=/Users/linmeixi/Lab/fin_whale/scripts_analyses/fsc26/param_btsp/fsc_output_list.csv
    WORKSCRIPT=/Users/linmeixi/Lab/fin_whale/scripts_analyses/fsc26/param_btsp/fsc26_generate_bootpar_${SETTING}.py
elif [ "${LOCAL}" == 'remote' ]; then
    set +xu
    source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
    conda activate gentools
    set -x
    SOURCE0=/u/project/rwayne/pnunez/Results/Demography/fsc2
    SINK0=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/fsc26/param_btsp/${SETTING}
    FSCFILE=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/fsc26/param_btsp/fsc_output_list.csv
    WORKSCRIPT=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/fsc26/param_btsp/fsc26_generate_bootpar_${SETTING}.py
else
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: Wrong LOCAL option ${SETTING}"
    exit 1
fi

# get the variables to use
MODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $2}' $FSCFILE)
SUBMODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $3}' $FSCFILE)
POP=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $4}' $FSCFILE)
PREFIX=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $5}' $FSCFILE)
DIR=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $6}' $FSCFILE)
EST=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $7}' $FSCFILE)
TPL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $8}' $FSCFILE)
PV=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $9}' $FSCFILE)
MAXL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $10}' $FSCFILE)

# the source and sink folders
SOURCE=${SOURCE0}/${DIR}
echo $SOURCE
ls $SOURCE
# check that source folder exists
if [ ! -d ${SOURCE} ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: No ${SOURCE}"
    exit 1
fi

SINK=${SINK0}/${MODEL}/${POP}
mkdir -p ${SINK}

# get the nchr variable for each populations
if [ "${POP}" == 'ENP' ]; then
    NCHR='3927079'
elif [ "${POP}" == 'GOC' ]; then
    NCHR='3908444'
elif [ "${POP}" == 'GOCENP' ]; then
    NCHR='3864185'
else
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: Wrong ${POP} specification"
    exit 1
fi

###########################################################
## main
# generate *_boot.par files
python ${WORKSCRIPT} ${SOURCE}/${MAXL} ${SINK}/${SUBMODEL}.${POP}_${PREFIX}_boot.par ${NCHR}

# copy and rename the est, tpl and pv files
rsync -ahvv ${SOURCE}/${EST} ${SINK}/${SUBMODEL}.${POP}_${PREFIX}_boot.est
rsync -ahvv ${SOURCE}/${TPL} ${SINK}/${SUBMODEL}.${POP}_${PREFIX}_boot.tpl
rsync -ahvv ${SOURCE}/${PV} ${SINK}/${SUBMODEL}.${POP}_${PREFIX}_boot.pv

echo -e "[$(date "+%Y-%m-%d %T")] Done"

