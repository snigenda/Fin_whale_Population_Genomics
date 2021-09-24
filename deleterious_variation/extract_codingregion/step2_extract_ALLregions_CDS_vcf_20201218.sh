#!/bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/gtf_CDS_step2_extract_ALLregions_CDS_vcf_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/gtf_CDS_step2_extract_ALLregions_CDS_vcf_20201218.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub step2_extract_ALLregions_CDS_vcf_20201218.sh
# @description	extract VCF from the CDS bedfiles
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 00:06:13 2020
# Run after: step1_output_CDS_from_gtf_20201217.sh

###########################################################
## import packages
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
# dataset
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")
# number of contigs
NJOBS=96

# directories
HOMEDIR="/u/project/rwayne/meixilin/fin_whale/analyses"
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

VCFDIR="/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}"
OUTDIR=${HOMEDIR}/get_ALLregions_CDS/${DATASET}/${REF}
LOGDIR=${OUTDIR}/logs
SCRATCHDIR="/u/scratch/m/meixilin/finwhale/analyses/get_ALLregions_CDS/${DATASET}/${REF}"

VCFPREFIX="JointCalls_${DATASET}_08_B_VariantFiltration"
OUTPREFIX="${VCFPREFIX}_${CDSTYPE}" # includes all the sites
OUTFILES=${LOGDIR}/${OUTPREFIX}_vcflist_${TODAY}.txt

# The bed file used
BEDFILE="/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.bed"
# The reference fasta
REFERENCE="/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta"

###########################################################
## main

### Get bed files and subset region for the selected type of variants
mkdir -p ${LOGDIR}
mkdir -p ${SCRATCHDIR}
cd ${SCRATCHDIR}
mkdir -p ./temp

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Subsetting vcf in ${VCFDIR} for Minke ALLregions CDS regions using ${BEDFILE} ..."

LOG=${LOGDIR}/gtf_CDS_step2_extract_ALLregions_CDS_vcf_${TODAY}.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Subsetting vcf in ${VCFDIR} for Minke ALLregions CDS regions using ${BEDFILE} ..." >> ${LOG}

for ii in $( seq 1 $NJOBS ); do
IDX=$(printf %02d ${ii})
MYVCF=${VCFDIR}/${VCFPREFIX}_${IDX}.vcf.gz
MYOUT=${SCRATCHDIR}/${OUTPREFIX}_${IDX}.vcf.gz
echo -e "[$(date "+%Y-%m-%d %T")] Working on ${MYVCF}"
echo -e "[$(date "+%Y-%m-%d %T")] Working on ${MYVCF}" >> ${LOG}

# IMPORTANT: use gatk SelectVariant over bcftools for increased robustness
# IMPORTANT: keep the original WGSproc8 output (some sites could be monomorphic while marked as 'VariantType=SNP' e.g. NW_006724458.1	218327)

gatk3 -Xmx8G -Djava.io.tmpdir=./temp \
-R ${REFERENCE} \
-T SelectVariants \
-V ${MYVCF} \
-o ${MYOUT} \
--preserveAlleles \
-L ${BEDFILE} &>> ${LOG}

echo ${MYOUT} >> ${OUTFILES}
echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
done

### Concatenate vcf files for SFS analyses
echo "[$(date "+%Y-%m-%d %T")] Concatenating vcf to ${OUTDIR}/${OUTPREFIX}_all.vcf.gz ..."
echo "[$(date "+%Y-%m-%d %T")] Concatenating vcf to ${OUTDIR}/${OUTPREFIX}_all.vcf.gz ..." >> ${LOG}

bcftools concat -f ${OUTFILES} -O z -o ${OUTDIR}/${OUTPREFIX}_all.vcf.gz
tabix -p vcf ${OUTDIR}/${OUTPREFIX}_all.vcf.gz

# generate the count
bcftools +counts ${OUTDIR}/${OUTPREFIX}_all.vcf.gz

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}

conda deactivate