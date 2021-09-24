#!/bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step3_snpEff_impact_extract_filteredvcf_vcf_20210129.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step3_snpEff_impact_extract_filteredvcf_vcf_20210129.err.txt
#$ -m abe

# @version        v0
# @usage          qsub step3_snpEff_impact_extract_filteredvcf_vcf_20210129.sh <impact_type>
# @description    Extract snpEff impact categories from filteredvcf (Note: --excludeFiltered --removeUnusedAlternates --excludeNonVariants)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Jan 31 15:20:31 2021
# Adapted from: step3_snpEff_LOF_extract_ALLregions_vcf_20201218.sh

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

extract_multiPASSSNPs_bedfile() {
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

# only allowing for PASS sites (should be already covered in the bedfile)
gatk3 -Xmx8G -Djava.io.tmpdir=./temp \
-R ${REFERENCE} \
-T SelectVariants \
--excludeFiltered \
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

echo -e "[$(date "+%Y-%m-%d %T")] Counting ${WORKDIR}/${OUTVCF}"
bcftools +counts ${WORKDIR}/${OUTVCF}

echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
}

###########################################################
## input variables
TYPE=${1} # MODIFIER/LOW/MODERATE/HIGH
DATASET="all50"
REF="Minke"
CDSTYPE="filteredvcf"
TODAY=$(date "+%Y%m%d")
NJOBS="96"

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
# directories
WORKDIR=${HOMEDIR}/snpEff_impact/${DATASET}/${REF}
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/snpEff_impact/${DATASET}/${REF}


# files
MYOUT=JointCalls_${DATASET}_filterpass_${TYPE}_${CDSTYPE}_snpEff
LOG=${WORKDIR}/logs/snpEff_${TYPE}_step3_extract_${CDSTYPE}_vcf_${TODAY}.log

# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

# input vcf
VCFDIR=/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}
VCFPREFIX=JointCalls_${DATASET}_08_B_VariantFiltration

###########################################################
## main
cd ${SCRATCHDIR}
mkdir -p ./temp

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}; Extracting ${MYOUT} from ${CDSTYPE}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}; Extracting ${MYOUT} from ${CDSTYPE}" > ${LOG}

extract_multiPASSSNPs_bedfile ${MYOUT}.vcf.gz ${WORKDIR}/bedfiles/${MYOUT}.bed

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}

conda deactivate
