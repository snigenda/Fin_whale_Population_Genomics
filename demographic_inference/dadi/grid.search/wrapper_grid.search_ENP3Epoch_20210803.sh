#!/bin/bash
#
# @version 		v1
# @script		bash wrapper_grid.search_ENP3Epoch_20210721.sh
# @description	local wrapper to run ENP 3Epoch grid search
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Jul 21 00:47:05 2021
# @modification: Tue Aug  3 15:59:29 2021
# @modification: Update to include Nanc values and more grid points

###########################################################
## import packages
source ~/miniconda2/etc/profile.d/conda.sh
conda activate py27
set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
# input scripts
TODAY=$(date "+%Y%m%d")
POP='ENP'
MODEL='1D.3Epoch'
SFS='/Users/linmeixi/Lab/fin_whale/scripts_analyses/dadi/grid.search/SFS/ENP-44.sfs'
L='392707916'
TBplusTF_Fix='0.134235802901809'
nuB_Fix='1.45113705136459'
nuF_Low='0.0001'
nuF_High='2'
TF_Low='1e-05'
TF_High='1'
numGridPoints='100'

# other
HOMEDIR=/Users/linmeixi/Lab/fin_whale/scripts_analyses
OUTDIR=/Users/linmeixi/google_drive/finwhale/analyses/dadi/grid.search/${POP}_${MODEL}_${TODAY}
LOG=${OUTDIR}/${POP}_${MODEL}_${TODAY}.log
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse master)
WORKSCRIPT=${HOMEDIR}/dadi/grid.search/grid.Search.1D.3Epoch.dadi.dadiUnits.py

###########################################################
## main
mkdir -p ${OUTDIR}
echo -e "[$(date "+%Y-%m-%d %T")] Running locally; GIT commit id ${COMMITID}" > ${LOG}

python ${WORKSCRIPT} \
--pop ${POP} \
--sfs $SFS \
--L ${L} \
--TBplusTF_Fix ${TBplusTF_Fix} \
--nuB_Fix ${nuB_Fix} \
--nuF_Low ${nuF_Low} \
--nuF_High ${nuF_High} \
--TF_Low ${TF_Low} \
--TF_High ${TF_High} \
--numGridPoints ${numGridPoints} \
--outdir ${OUTDIR} >> ${LOG} 2>&1
exitVal=${?}
echo "python exit Value = ${exitVal}"
echo "python exit Value = ${exitVal}" >> ${LOG}

conda list >> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] Everything Done" >> ${LOG}

conda deactivate
