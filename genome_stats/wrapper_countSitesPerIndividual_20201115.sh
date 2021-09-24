#!/bin/bash
#$ -l h_data=4G,h_vmem=10G,h_rt=20:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -N count_sites_20201115
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/count_sites_20201115.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/count_sites_20201115.err.txt
#$ -m abe

# @version 		v1
# @script		countSitesPerIndividual_record.sh
# @description	Wrapper of calling the countSitesPerIndividual.py in Hoffman2
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon May 11 00:51:52 2020
# @modification Fri Aug 28 12:01:49 2020
# @modification Update to match all50 dataset
# @modification: Sun Nov 15 15:37:33 2020
# @modification: Update to match Nov 2020 version of all50 dataset

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
DATASET="all50"
REF="Minke"

## def variables
TODAY=$(date "+%Y%m%d")
WORK=/u/project/rwayne/meixilin/fin_whale/analyses
WORKSCRIPT=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/Summary_stats/count_sites/countSitesPerIndividual.py
OUTDIR=${WORK}/Summary_stats/${DATASET}/${REF}/count_sites_${TODAY}
VCFDIR=/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}

IDX=$(printf %02d ${SGE_TASK_ID})
COMMITID=$(git --git-dir="${WORK}/scripts/.git" --work-tree="${WORK}/scripts" rev-parse master)
VCF="JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz"

mkdir -p ${OUTDIR}

###########################################################
## main
# first check the overall vcf.gz file
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Count sites for ${VCF}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

python ${WORKSCRIPT} --vcf ${VCFDIR}/${VCF} --outfile ${OUTDIR}/${DATASET}_${REF}_${IDX}_sites_summary.txt --filter "PASS" --contig ${IDX}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"

conda deactivate
