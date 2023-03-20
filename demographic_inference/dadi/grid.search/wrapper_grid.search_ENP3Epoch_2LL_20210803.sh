#!/bin/bash
#
# @version 		v1
# @script		bash wrapper_grid.search_ENP3Epoch_1LL_20210803.sh
# @description	local wrapper to run ENP 3Epoch grid search. Zoom in on focal regions
# Author: Meixi Lin
# Date: Wed Jul 21 00:47:05 2021
# @modification: Wed Aug  4 11:36:58 2021
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
SFS='<homedir>/fin_whale/scripts_analyses/dadi/grid.search/SFS/ENP-44.sfs'
L='392707916'
TBplusTF_Fix='0.134235802901809'
nuB_Fix='1.45113705136459'
nuF_Low='0.0030000010909600'
nuF_High='0.3651483052764145'
TF_Low='0.0000100000000000'
TF_High='0.0016213206814270'
intercept_Low='-2.604425'
intercept_High='-2.204425'
slope='1.040538'
numGridPoints='400'

# other
HOMEDIR=<homedir>/fin_whale/scripts_analyses
OUTDIR=<homedir>/finwhale/analyses/dadi/grid.search/${POP}_${MODEL}_2LL_${TODAY}
LOG=${OUTDIR}/${POP}_${MODEL}_2LL_${TODAY}.log
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse master)
WORKSCRIPT=${HOMEDIR}/dadi/grid.search/grid.Search.Ridge.1D.3Epoch.dadi.dadiUnits.py

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
--intercept_Low ${intercept_Low} \
--intercept_High ${intercept_High} \
--slope ${slope} \
--numGridPoints ${numGridPoints} \
--outdir ${OUTDIR} >> ${LOG} 2>&1
exitVal=${?}
echo "python exit Value = ${exitVal}"
echo "python exit Value = ${exitVal}" >> ${LOG}

conda list >> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] Everything Done" >> ${LOG}

conda deactivate
