#!/bin/bash
#$ -l h_data=4G,h_vmem=8G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_nonsyn_type_step1_extract_ALLregions_bed_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_nonsyn_type_step1_extract_ALLregions_bed_20201218.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub step1_extract_NONSYNtype_ALLregions_bed_20201218.sh
# @description	extract DELETERIOUS and TOLERATED annotations
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Dec 19 14:32:28 2020
# @modification: Mon Dec 28 13:54:49 2020
# @modification: New nonsyn regions with triming alternatives applied

###########################################################
## import packages
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")
PASSINGVAR="PASS,WARN_missing"

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
WORKSCRIPT=${HOMEDIR}/scripts/get_ALLregions_CDS/SIFT_nonsyn_type/extract_NONSYNtype_SIFT_bed.py

# directories
WORKDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}
OUTDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles

# input vcf
MYVCF=${WORKDIR}/JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.vcf.gz
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Extracting NONSYN type from ALregions, SIFT in ${MYVCF}"

cd ${WORKDIR}

echo -e "COMMAND: python $WORKSCRIPT --VCF ${MYVCF} --filter ${PASSINGVAR} --outprefix ${MYVCF/.vcf.gz}"
python $WORKSCRIPT --VCF ${MYVCF} --filter ${PASSINGVAR} --outprefix ${MYVCF/.vcf.gz}

mv -v ${MYVCF/.vcf.gz/_NONSYNtype.bed} ${OUTDIR}/

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
conda deactivate