#!/bin/bash
#
# @version 		v0
# @script		bash step3_Prunedvcf2plink_all50_20210316.sh
# @description	convert a given vcf file to plink format
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Mar 18 14:42:13 2021

###########################################################
## import packages
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -xeo pipefail

###########################################################
## def functions

###########################################################
## def variables
VCFFILE=${1}
OUTPREFIX=${2}
###########################################################
## main
# The '12' modifier causes A1 (usually minor) alleles to be coded as '1'
#       and A2 alleles to be coded as '2', while '01' maps A1 -> 0 and A2 -> 1.
plink --vcf ${VCFFILE} --allow-extra-chr --recode 12 --out ${OUTPREFIX}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
