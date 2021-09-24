#!/bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_LOF_step1_extract_ALLregions_bed_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_LOF_step1_extract_ALLregions_bed_20201218.err.txt
#$ -m abe

# @version        v0
# @usage          qsub -t 1-96 step1_snpEff_LOF_extract_ALLregions_bed_20201218.sh
# @description    Add loss-of-function annotations bedfiles using snpEff
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 13:20:48 2020
# Criteria: http://snpeff.sourceforge.net/snpEff_lof_nmd.pdfs

###########################################################
## import packages
sleep $((RANDOM % 120))

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
IDX=$(printf %02d ${SGE_TASK_ID})
TODAY=$(date "+%Y%m%d")

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKSCRIPT=${HOMEDIR}/scripts/get_ALLregions_CDS/snpEff_LOF/extract_LOFregion_snpEff_bed.py
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
VCFDIR="/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}"
OUTDIR=${HOMEDIR}/get_ALLregions_CDS//${DATASET}/${REF}/bedfiles
LOGDIR=${HOMEDIR}/get_ALLregions_CDS//${DATASET}/${REF}/logs
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles
mkdir -p ${OUTDIR}
mkdir -p ${SCRATCHDIR}


# input vcf
MYVCF=${VCFDIR}/JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## main

cd ${SCRATCHDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

PASSINGVAR="PASS,WARN_missing"
CUTOFF="0.5"

# vcf + dataset + filter + LOF cutoff + cdstype + chromosome section + software
MYOUT=JointCalls_${DATASET}_filterpassmiss_LOF05_${CDSTYPE}_${IDX}_snpEff

echo -e "[$(date "+%Y-%m-%d %T")] Extracting LOF from ${MVCF} to ${SCRATCHDIR}/${MYOUT}.bed, cutoff = ${CUTOFF}, Passing sites = ${PASSINGVAR}"

python $WORKSCRIPT --VCF ${MYVCF} --filter ${PASSINGVAR} --cutoff ${CUTOFF} --outprefix ${MYOUT}
check_length ${MYOUT}.bed
echo -e "[$(date "+%Y-%m-%d %T")] Done"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"

conda deactivate
