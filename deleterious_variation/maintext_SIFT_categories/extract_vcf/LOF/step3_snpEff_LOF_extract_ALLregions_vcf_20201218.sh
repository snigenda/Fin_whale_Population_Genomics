#!/bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_LOF_step3_extract_ALLregions_vcf_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/snpEff_LOF_step3_extract_ALLregions_vcf_20201218.err.txt
#$ -m abe

# @version        v1
# @usage          qsub step3_snpEff_LOF_extract_ALLregions_vcf_20201218.sh
# @description    Extract LOF vcf files from ALLregions
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Dec  5 16:20:11 2020
# @modification: Mon Dec 28 15:30:45 2020
# @modification: New LOF regions with triming alternatives applied

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

extract_multiSNPs_bedfile() {
local OUTVCF=${1}
local BEDFILE=${2}

local OUTVCFLIST=${SCRATCHDIR}/${OUTVCF/.vcf.gz/_vcflist.txt}
echo -e "[$(date "+%Y-%m-%d %T")] Extracting vcf from ${BEDFILE} outputting to ${WORKDIR}/${OUTVCF}"
echo -e "[$(date "+%Y-%m-%d %T")] Extracting vcf from ${BEDFILE} outputting to ${WORKDIR}/${OUTVCF}" >> ${LOG}

# loop through the vcfs
for ii in $( seq 1 $NJOBS ); do
local IDX=$(printf %02d ${ii})
local MYVCF=${VCFDIR}/${VCFPREFIX}_${IDX}.vcf.gz
local OUTVCFIDX=${SCRATCHDIR}/${OUTVCF/.vcf.gz/_${IDX}.vcf.gz}

echo -e "[$(date "+%Y-%m-%d %T")] Working on ${MYVCF} outputting to ${OUTVCFIDX}" >> ${LOG}

gatk3 -Xmx8G -Djava.io.tmpdir=./temp \
-R ${REFERENCE} \
-T SelectVariants \
--removeUnusedAlternates \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-V ${MYVCF} \
-o ${OUTVCFIDX} \
-L ${BEDFILE} &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi
echo ${OUTVCFIDX} >> ${OUTVCFLIST}
done

bcftools concat -f ${OUTVCFLIST} -O z -o ${WORKDIR}/${OUTVCF}
tabix -p vcf ${WORKDIR}/${OUTVCF}

bcftools +counts ${WORKDIR}/${OUTVCF}

echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
}

###########################################################
## input variables
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")
NJOBS="96"

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/get_ALLregions_CDS/${DATASET}/${REF}
WORKDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}

BEDDIR=${WORKDIR}/bedfiles
LOGDIR=${WORKDIR}/logs
LOG=${LOGDIR}/snpEff_LOF_step3_extract_ALLregions_vcf_${TODAY}.log

# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

# input vcf
VCFDIR=/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}
VCFPREFIX="JointCalls_${DATASET}_08_B_VariantFiltration"

###########################################################
## main
cd ${SCRATCHDIR}
mkdir -p ./temp

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${LOG}

MYOUT=JointCalls_${DATASET}_filterpassmiss_LOF05_${CDSTYPE}_all_snpEff
extract_multiSNPs_bedfile ${MYOUT}.vcf.gz ${BEDDIR}/${MYOUT}.bed

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}


conda deactivate
