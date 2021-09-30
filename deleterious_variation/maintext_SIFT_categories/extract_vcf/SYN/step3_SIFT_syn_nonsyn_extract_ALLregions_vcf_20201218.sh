#!/bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_syn_nonsyn_step3_extract_ALLregions_vcf_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_syn_nonsyn_step3_extract_ALLregions_vcf_20201218.err.txt
#$ -m abe

# @version      v2
# @usage        qsub step3_SIFT_syn_nonsyn_extract_ALLregions_vcf_20201218.sh
# @description  Wrapper of extracting for SYN/NONSYN mutations in ALLregions (identified by SIFT) from bedfiles as vcf.gz in Hoffman2
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 20:43:57 2020
# should be executed after: step2_SIFT_syn_nonsyn_merge_ALLregions_bed_20201218.sh
# @modification: Fri Dec 18 20:57:52 2020
# @modification: Update to exclude sites marked as SNPs in the bedfiles but were monomorphic (e.g. the variant sites did not pass the genotype filter)

###########################################################
## import packages
sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail # for safer scripting

###########################################################
## def functions
# extract vcf files according to bedfiles
extract_SNPs_bedfile() {
local MYVCF=${1}
local OUTVCF=${2}
local BEDFILE=${3}

echo -e "[$(date "+%Y-%m-%d %T")] Extracting ${MYVCF} from ${BEDFILE} outputting to ${OUTVCF}"

gatk3 -Xmx8G -Djava.io.tmpdir=./temp \
-R ${REFERENCE} \
-T SelectVariants \
--removeUnusedAlternates \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-V ${MYVCF} \
-o ${OUTVCF} \
-L ${BEDFILE}

bcftools +counts ${OUTVCF}

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
OUTDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}
BEDDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles

mkdir -p ${OUTDIR}

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

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] outputting SYN/NONSYN vcf subset from ${MYVCF}; PassingVar = ${PASSINGVAR}; removeUnusedAlternates."

PASSINGVAR="PASS,WARN_missing"
# for SIFT
SYN_SIFT="SYNONYMOUS"
NONSYN_SIFT="NONSYNONYMOUS"

###########################################################
# 1. SYN
MYOUT=JointCalls_all50_filterpassmiss_syn_${CDSTYPE}_all_SIFT
extract_SNPs_bedfile ${MYVCF} ${MYOUT}.vcf.gz ${BEDDIR}/${MYOUT}.bed

####################e#######################################
# 2. NONSYN
MYOUT=JointCalls_all50_filterpassmiss_nonsyn_${CDSTYPE}_all_SIFT
extract_SNPs_bedfile ${MYVCF} ${MYOUT}.vcf.gz ${BEDDIR}/${MYOUT}.bed

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate
