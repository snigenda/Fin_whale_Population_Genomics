#!/bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step1_snpEff_impact_extract_filteredvcf_bed_20210129.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step1_snpEff_impact_extract_filteredvcf_bed_20210129.err.txt
#$ -m abe

# @version        v0
# @usage          qsub -t 1-96 step1_snpEff_impact_extract_filteredvcf_bed_20210129.sh
# @description    Traverse and extract the Annotation_Impact from snpEff Annotation
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jan 29 15:18:17 2021

###########################################################
## import packages
sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail # for safer scripting

###########################################################
## def functions

###########################################################
## input variables
DATASET="all50"
REF="Minke"
CDSTYPE="filteredvcf"
IDX=$(printf %02d ${SGE_TASK_ID})
TODAY=$(date "+%Y%m%d")

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKSCRIPT=${HOMEDIR}/scripts/snpEff_impact/extract_Annotation_Impact_snpEff_bed.py
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
VCFDIR=/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/snpEff_impact/${DATASET}/${REF}/bedfiles
mkdir -p ${SCRATCHDIR}

# input vcf
MYVCF=${VCFDIR}/JointCalls_${DATASET}_08_B_VariantFiltration_${IDX}.vcf.gz
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

# variables
PASSINGVAR="PASS"
# vcf + dataset + filter + analyses type + cdstype + chromosome section + software
MYOUT=JointCalls_${DATASET}_filterpass_impact_${CDSTYPE}_${IDX}_snpEff

###########################################################
## mains
cd ${SCRATCHDIR}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] COMMAND: python $WORKSCRIPT --VCF ${MYVCF} --filter ${PASSINGVAR} --outprefix ${MYOUT}"

python $WORKSCRIPT --VCF ${MYVCF} --filter ${PASSINGVAR} --outprefix ${MYOUT}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"

conda deactivate
