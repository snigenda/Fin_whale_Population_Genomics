#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/Summary_stats/f50b4_countSitesPerIndividual_20210127.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/Summary_stats/f50b4_countSitesPerIndividual_20210127.err.txt
#$ -m abe

# @version      v1
# @script       wrapper_f50b4_countSitesPerIndividual_20210127.sh
# @usage        qsub -t 1-96 wrapper_f50b4_countSitesPerIndividual_20210127.sh
# @description  Wrapper of calling the countSitesPerIndividual.py for the baleen_genomes dataset
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Jan 27 13:24:01 2021
# Adapted from wrapper_ALLregions_CDS_countSitesPerIndividual_20201228.sh

###########################################################
## import packages
sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools
set -eo pipefail

###########################################################
## def functions

###########################################################
## input variables
DATASET="f50b4"
REF="Minke"
CDSTYPE="filteredvcf"
IDX=$(printf %02d ${SGE_TASK_ID})

## def variables
TODAY=$(date "+%Y%m%d")
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKSCRIPT=${HOMEDIR}/scripts/Summary_stats/count_sites/countSitesPerIndividual.py
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

VCFDIR=${HOMEDIR}/baleen_genomes/filteredvcf/${DATASET}/${REF}
OUTDIR=${HOMEDIR}/Summary_stats/${DATASET}/${REF}/count_sites_${CDSTYPE}_${TODAY}
VCF="JointCalls_${DATASET}_08_B_VariantFiltration_${IDX}.vcf.gz"

mkdir -p ${OUTDIR}

###########################################################
## main
# first check the overall vcf.gz file
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; Count sites for ${VCFDIR}/${VCF}; WORKSCRIPT=${WORKSCRIPT}; git commit id: ${COMMITID}"

python ${WORKSCRIPT} --vcf ${VCFDIR}/${VCF} --outfile ${OUTDIR}/${DATASET}_${REF}_${IDX}_sites_summary.txt --filter "PASS" --contig ${IDX}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"

conda deactivate
