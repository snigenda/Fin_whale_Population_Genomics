#!/bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -l h_data=4G,h_vmem=8G,h_rt=23:00:00
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_LOF_step2_merge_ALLregions_bed_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_LOF_step2_merge_ALLregions_bed_20201218.err.txt
#$ -m abe

# @version        v0
# @usage          qsub step2_snpEff_LOF_merge_ALLregions_bed_20201218.sh
# @description    Merge loss-of-function annotations bedfiles using snpEff
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Dec 28 15:20:44 2020
# Should be run after step1_snpEff_LOF_extract_ALLregions_bed_20201218.sh

###########################################################
## import packages

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail # for safer scripting

###########################################################
## def functions
check_length() {
    local FILENAME=${1}
    echo -e "[$(date "+%Y-%m-%d %T")] Getting total length of ${FILENAME}"
    awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${FILENAME}
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
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles
OUTDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles

# output prefix
MYOUT=JointCalls_${DATASET}_filterpassmiss_LOF05_${CDSTYPE}_all_snpEff
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## main
cd ${OUTDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Merging ${OUTDIR}/${MYOUT}.bed from ${SCRATCHDIR}/JointCalls_all50_filterpassmiss_LOF05_ALLregions_01-96_snpEff.bed"

cat ${SCRATCHDIR}/JointCalls_${DATASET}_filterpassmiss_LOF05_${CDSTYPE}_[0-9][0-9]_snpEff.bed | \
bedtools sort -i stdin > ${MYOUT}.bed

check_length ${MYOUT}.bed

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate
