#!/bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_nonsyn_type_step2_extract_ALLregions_vcf_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_nonsyn_type_step2_extract_ALLregions_vcf_20201218.err.txt
#$ -m abe

# @version      v1
# @usage        qsub step2_extract_NONSYNtype_ALLregions_vcf_20201218.sh
# @description  extract NONSYNtype from JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.vcf.gz
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Dec 14 23:08:34 2020
# Should be run after: step1_extract_NONSYNtype_ALLregions_bed_20201218.sh
# @modification: Mon Dec 28 14:43:40 2020
# @modification: New nonsyn regions with triming alternatives applied

###########################################################
## import packages

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions
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
## def variables
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")
PASSINGVAR="PASS,WARN_missing"

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
WORKDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}
BEDDIR=${WORKDIR}/bedfiles

# input vcf
MYVCF=${WORKDIR}/JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.vcf.gz
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Extracting ${MYVCF} by NONSYN types (TOLERATED/DELETERIOUS)"

cd ${WORKDIR}
mkdir -p ./temp
# DELETERIOUS with no warning present
MYOUT=JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_DELnoW
extract_SNPs_bedfile ${MYVCF} ${MYOUT}.vcf.gz ${BEDDIR}/${MYOUT}.bed

# DELETERIOUS (can have warning)
MYOUT=JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_DEL
extract_SNPs_bedfile ${MYVCF} ${MYOUT}.vcf.gz ${BEDDIR}/${MYOUT}.bed

# TOLERATED
MYOUT=JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_TOL
extract_SNPs_bedfile ${MYVCF} ${MYOUT}.vcf.gz ${BEDDIR}/${MYOUT}.bed

rm -r ./temp
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate