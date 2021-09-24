#!/bin/bash
#$ -l h_data=20G,h_vmem=24G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_impactCDS_step1_extract_ALLregions_vcf_20210415.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_impactCDS_step1_extract_ALLregions_vcf_20210415.err.txt
#$ -m abe

# @version      v0
# @usage        qsub step1_snpEff_impactCDS_extract_ALLregions_vcf_20210415.sh
# @description  Wrapper of extracting for HIGH/MODERATE/LOW/MODIFIER mutations in ALLregions (identified by snpEff) from bedfiles as vcf.gz in Hoffman2
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Apr 15 15:56:27 2021
# Note: The bedfiles were generated from: snpEff_impact/step1_snpEff_impact_extract_filteredvcf_bed_20210129.sh and snpEff_impact/step2_snpEff_impact_merge_filteredvcf_bed_20210129.sh (The bedfile included all the sites that had genotypes called)
# Additional filters: 1. excludeFiltered 2. excludeNonVariants

###########################################################
## import packages
sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -xeo pipefail # for safer scripting

###########################################################
## def functions
# extract vcf files according to bedfiles
extract_SNPs_bedfile() {
local MYVCF=${1}
local OUTVCF=${2}
local BEDFILE=${3}

echo -e "[$(date "+%Y-%m-%d %T")] Extracting ${MYVCF} from ${BEDFILE} outputting to ${OUTVCF}"

gatk3 -Xmx15G -Djava.io.tmpdir=./temp \
-R ${REFERENCE} \
-T SelectVariants \
--excludeFiltered \
--removeUnusedAlternates \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-V ${MYVCF} \
-o ${OUTVCF} \
-L ${BEDFILE} &>> ${LOG}

bcftools +counts ${OUTVCF}

echo -e "[$(date "+%Y-%m-%d %T")] Done"
}

###########################################################
## input variables
TYPE=${1} # MODIFIER/LOW/MODERATE/HIGH
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
PASSINGVAR="PASS"
TODAY=$(date "+%Y%m%d")

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
OUTDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}
BEDDIR=${HOMEDIR}/snpEff_impact/${DATASET}/${REF}/bedfiles
LOG=${OUTDIR}/logs/snpEff_impactCDS_${TYPE}_step1_extract_ALLregions_vcf_${TODAY}.log

# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

# input vcf
MYVCF=${OUTDIR}/JointCalls_all50_08_B_VariantFiltration_ALLregions_all.vcf.gz
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## subset vcfs

cd ${OUTDIR}
mkdir -p ./temp

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}; Type = ${TYPE}; PassingVar = ${PASSINGVAR}; removeUnusedAlternates"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}; Type = ${TYPE}; PassingVar = ${PASSINGVAR}; removeUnusedAlternates" > ${LOG}


MYOUT=JointCalls_${DATASET}_filterpass_${TYPE}_${CDSTYPE}_all_snpEff
MYBED=JointCalls_${DATASET}_filterpass_${TYPE}_filteredvcf_snpEff

extract_SNPs_bedfile ${MYVCF} ${OUTDIR}/${MYOUT}.vcf.gz ${BEDDIR}/${MYBED}.bed

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}

conda deactivate
