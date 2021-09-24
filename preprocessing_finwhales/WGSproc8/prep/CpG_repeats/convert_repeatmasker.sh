#!/bin/bash
#$ -l h_data=4G,h_rt=23:00:00
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -o /u/project/rwayne/snigenda/finwhale/reports/convert_repeatmasker_20200728.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/convert_repeatmasker_20200728.err.txt
#$ -m abe

# 
# @version      v0
# @script       convert_repeatmasker.sh
# @description  conver the repeat masker output to bed file (1-based to 0-based)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Jul 28 00:12:08 2020

###########################################################
## import packages 
set -o pipefail 

###########################################################
## def functions 

###########################################################
## def variables 
RMOUT="/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_rm.out"
RMOUT2=${RMOUT}.v2.bed

###########################################################
## main 

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"  
cd /u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0

awk  -F ' ' 'BEGIN{ OFS = "\t" }{
    if (NR < 4) {}
    else {$6 = $6 - 1 ; print $5,$6,$7} 
}' ${RMOUT} > ${RMOUT2}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" 
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" 