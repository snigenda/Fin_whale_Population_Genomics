#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=20:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_syn_nonsyn_step2_merge_ALLregions_bed_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_syn_nonsyn_step2_merge_ALLregions_bed_20201218.err.txt
#$ -m abe

# @version      v1
# @script       qsub step2_SIFT_syn_nonsyn_merge_ALLregions_bed_20201218.sh
# @description  Mergeing the output of step1_SIFT_syn_nonsyn_extract_ALLregions_bed_20201218.sh
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 11:13:12 2020

###########################################################
## import packages
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail # for safer scripting

###########################################################
## def functions
# check length of file
check_length() {
    local FILENAME=${1}
    echo -e "[$(date "+%Y-%m-%d %T")] Getting total length of ${FILENAME}"
    awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${FILENAME}
}

# write length of file
merge_bedfiles() {
    local PREFIX=${1}
    local ANNTYPE=${2}
    local SOFTWARE=${3}

    echo -e "[$(date "+%Y-%m-%d %T")] Merging ${OUTDIR}/${PREFIX}.bed from ${SCRATCHDIR}/JointCalls_all50_filterpassmiss_${ANNTYPE}_ALL_01-96_${SOFTWARE}.bed"

    cat ${SCRATCHDIR}/JointCalls_all50_filterpassmiss_${ANNTYPE}_ALLregions_[0-9][0-9]_${SOFTWARE}.bed | \
    bedtools sort -i stdin > ${OUTDIR}/${PREFIX}.bed
    check_length ${OUTDIR}/${PREFIX}.bed

    echo -e "[$(date "+%Y-%m-%d %T")] Done"
}

###########################################################
## input variables
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
OUTDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles
LOGDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}/logs
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles

mkdir -p ${SCRATCHDIR}
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## logging

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Merging SYN/NONSYN from ALLregions CDS in ${SCRATCHDIR}"

###########################################################
## variables used in the extraction before
cd ${OUTDIR}

PASSINGVAR="PASS,WARN_missing"

# for SIFT
SYN_SIFT="SYNONYMOUS"
NONSYN_SIFT="NONSYNONYMOUS"

###########################################################
# 1. SYN
# Merge the bedfiles
MYANN="syn"
MYSOFT="SIFT"
MYOUT=JointCalls_all50_filterpassmiss_${MYANN}_${CDSTYPE}_all_${MYSOFT}
merge_bedfiles ${MYOUT} ${MYANN} ${MYSOFT}

###########################################################
# 2. NONSYN
MYANN="nonsyn"

MYSOFT="SIFT"
MYOUT=JointCalls_all50_filterpassmiss_${MYANN}_${CDSTYPE}_all_${MYSOFT}
merge_bedfiles ${MYOUT} ${MYANN} ${MYSOFT}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate
