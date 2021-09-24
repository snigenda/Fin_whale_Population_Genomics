#!/bin/bash
#$ -l h_data=20G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step4_admixture_all50_20210316.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step4_admixture_all50_20210316.err.txt
#$ -m abe

# @version      v0
# @usage        qsub step4_admixture_all50_20210316.sh
# @description  Performs admixture analyses on the LDpruned sites (MAF cutoff = NA/0.05/0.10)
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Thu Mar 18 15:22:26 2021
# Output:
# 1) 10 runs of admixture from K=2 to K=6
# 2) files with stats
# Notes:
# inputFile may be:
#      - a PLINK .bed file
#      - a PLINK "12" coded .ped file

############################################################
## import packages
set -eo pipefail

############################################################
## def functions

############################################################
## def variables
DATASET='all50'
REF='Minke'
MAFLISTS=('05' '10' 'NA')
TODAY=$(date "+%Y%m%d")

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}


COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# admixture software source: http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
ADMIXTURE=/u/project/rwayne/meixilin/software/admixture/dist/admixture_linux-1.3.0/admixture

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}"

for MAFCUT in ${MAFLISTS[@]}; do
    echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Start."
    OUTDIR=${WORKDIR}/Admixture_${TODAY}/maf${MAFCUT}
    mkdir -p ${OUTDIR}
    cd ${OUTDIR}
    # filenames
    OUT=JointCalls_${DATASET}_filterpass_bialleic_all_LDPruned_maf${MAFCUT}_SA_mrF
    SHORTOUT=${DATASET}_pass_maf${MAFCUT}
    PED=${WORKDIR}/${OUT}.ped
    echo -e "K,iter,CVERROR,LL" > Admixture_CV_LLsummary_maf${MAFCUT}.csv
    for K in {2..6};do
        for i in {1..10};do
            # -s time the random seed to be generated from the current time
            # --cv In this default setting, the cross-validation procedure will perform 5-fold CV
            ${ADMIXTURE} --cv -s time -j8 ${PED} ${K} | tee log_K${K}.iter${i}.out
            mv ${OUT}.${K}.Q ${SHORTOUT}.K${K}.iter${i}.Q
            mv ${OUT}.${K}.P ${SHORTOUT}.K${K}.iter${i}.P
            # get the CV error and loglikelihood during each run
            CVERROR=$(awk '/^CV/ {print $4}' log_K${K}.iter${i}.out)
            LL=$(awk '/^Loglikelihood/ {print $2}' log_K${K}.iter${i}.out)
            echo -e "${K},${i},${CVERROR},${LL}" >> Admixture_CV_LLsummary_maf${MAFCUT}.csv
        done
    done
    echo -e "[$(date "+%Y-%m-%d %T")] Maf cutoff = ${MAFCUT}; Done."
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
