#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=20:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_syn_nonsyn_step1_extract_ALLregions_bed_20201218.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/SIFT_syn_nonsyn_step1_extract_ALLregions_bed_20201218.err.txt
#$ -m abe

# @version      v1
# @usage        qsub -t 1-96 step1_SIFT_syn_nonsyn_extract_ALLregions_bed_20201218.sh
# @description  Wrapper of extracting for syn/nonsyn mutation regions from ALL the regions called as bedfiles in Hoffman2 (to match up the vaquita analyses) using the SIFT annotations
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 01:01:15 2020

###########################################################
## import packages
sleep $((RANDOM % 120))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail # for safer scripting

###########################################################
## def functions
# write length of file
check_length() {
    local FILENAME=${1}
    echo -e "[$(date "+%Y-%m-%d %T")] Getting total length of ${FILENAME}"
    awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${FILENAME}
}

# extract the given type of annotations
extract_annotation_bed_SIFT() {
    local ANNTYPE=${1}
    local PREFIX=${2}
    local VERBOSE=${3}

    echo -e "[$(date "+%Y-%m-%d %T")] Extracting ${SCRATCHDIR}/${PREFIX}.bed from ${MYVCF}, Annotation=${ANNTYPE}, Verbosity=${VERBOSE}"
    echo -e "[$(date "+%Y-%m-%d %T")] Extracting ${SCRATCHDIR}/${PREFIX}.bed from ${MYVCF}, Annotation=${ANNTYPE}, Verbosity=${VERBOSE}" >> ${LOG}

    if [[ ${VERBOSE} == "True" ]]; then
        echo -e "COMMAND: python $WORKSCRIPT_SIFT --VCF ${MYVCF} --anntype_SIFT ${ANNTYPE} --filter ${PASSINGVAR} --outprefix ${PREFIX} --verbose 2> ${LOGDIR}/${PREFIX}.verbose.err.${TODAY}.txt" >> ${LOG}
        python $WORKSCRIPT_SIFT --VCF ${MYVCF} --anntype_SIFT ${ANNTYPE} --filter ${PASSINGVAR} --outprefix ${PREFIX} --verbose 2> ${LOGDIR}/${PREFIX}.verbose.err.${TODAY}.txt
        gzip ${LOGDIR}/${PREFIX}.verbose.err.${TODAY}.txt
    else
        echo -e "COMMAND: python $WORKSCRIPT_SIFT --VCF ${MYVCF} --anntype_SIFT ${ANNTYPE} --filter ${PASSINGVAR} --outprefix ${PREFIX}" >> ${LOG}
        python $WORKSCRIPT_SIFT --VCF ${MYVCF} --anntype_SIFT ${ANNTYPE} --filter ${PASSINGVAR} --outprefix ${PREFIX} 2>> ${LOG}
    fi

    check_length ${PREFIX}.bed

    echo -e "[$(date "+%Y-%m-%d %T")] Done"
    echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
}

###########################################################
## input variables
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")
IDX=$(printf %02d ${SGE_TASK_ID})

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
WORKSCRIPT_SIFT=${HOMEDIR}/scripts/get_ALLregions_CDS/SIFT_syn_nonsyn/extract_annregion_SIFT_bed.py
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# directories
VCFDIR="/u/project/rwayne/snigenda/finwhale/filteredvcf/${DATASET}/${REF}"
OUTDIR=${HOMEDIR}/get_ALLregions_CDS//${DATASET}/${REF}/bedfiles
LOGDIR=${HOMEDIR}/get_ALLregions_CDS//${DATASET}/${REF}/logs
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/get_ALLregions_CDS/${DATASET}/${REF}/bedfiles

# input vcf
MYVCF=${VCFDIR}/JointCalls_${DATASET}_08_B_VariantFiltration_${IDX}.vcf.gz
# reference genome
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

###########################################################
## logging
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}
mkdir -p ${SCRATCHDIR}
cd ${SCRATCHDIR}
mkdir -p ./temp

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Extracting SYN/NONSYN from ALLregions CDS in ${MYVCF}"

LOG=${LOGDIR}/SIFT_syn_nonsyn_step1_extract_ALLregions_bed_chr${IDX}_${TODAY}.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Extracting SYN/NONSYN from ALLregions CDS in ${MYVCF}" >> ${LOG}

###########################################################
## extract bedfiles for different types of mutations
cd ${SCRATCHDIR}

PASSINGVAR="PASS,WARN_missing"

# for SIFT
SYN_SIFT="SYNONYMOUS"
NONSYN_SIFT="NONSYNONYMOUS"

###########################################################
# 1. SYN
# be verbose and include information about the transcript chosen
MYOUT=JointCalls_${DATASET}_filterpassmiss_syn_${CDSTYPE}_${IDX}_SIFT
extract_annotation_bed_SIFT ${SYN_SIFT} ${MYOUT} "True"

###########################################################
# 2. NONSYN
MYOUT=JointCalls_${DATASET}_filterpassmiss_nonsyn_${CDSTYPE}_${IDX}_SIFT
extract_annotation_bed_SIFT ${NONSYN_SIFT} ${MYOUT} "False"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}

conda deactivate
