#!/bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -l h_data=4G,h_vmem=8G,h_rt=23:00:00
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step2_snpEff_impact_merge_filteredvcf_bed_20210129.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step2_snpEff_impact_merge_filteredvcf_bed_20210129.err.txt
#$ -m abe

# @version        v0
# @usage          qsub step2_snpEff_impact_merge_filteredvcf_bed_20210129.sh
# @description    Merge snpEff Annotation_Impact bedfiles into three categories
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Dec 28 15:20:44 2020
# Should be run after step1_snpEff_LOF_extract_ALLregions_bed_20201218.sh

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
    if [[ ${FILENAME} == *.gz ]]; then
        zcat ${FILENAME} | awk -F'\t' '{SUM += $3-$2} END {print SUM}'
    elif [[ ${FILENAME} == *.bed ]]; then
       awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${FILENAME}
    fi
}

# write tally of length
output_tally() {
    local LINES=$(wc -l ${BEDFILE}|cut -d' ' -f 1)
    local COUNT1=$(awk 'BEGIN {FS = "\t"}; $4 == "MODIFIER" {count ++}; END {print count}' ${BEDFILE})
    local COUNT2=$(awk 'BEGIN {FS = "\t"}; $4 == "LOW" {count ++}; END {print count}' ${BEDFILE})
    local COUNT3=$(awk 'BEGIN {FS = "\t"}; $4 == "MODERATE" {count ++}; END {print count}' ${BEDFILE})
    local COUNT4=$(awk 'BEGIN {FS = "\t"}; $4 == "HIGH" {count ++}; END {print count}' ${BEDFILE})
    echo -e "${IDX}\t${LINES}\t${COUNT1}\t${COUNT2}\t${COUNT3}\t${COUNT4}"
    echo -e "${IDX}\t${LINES}\t${COUNT1}\t${COUNT2}\t${COUNT3}\t${COUNT4}" >> ${LOG}
}

###########################################################
## input variables
DATASET="all50"
REF="Minke"
CDSTYPE="filteredvcf"
NJOBS="96"
TODAY=$(date "+%Y%m%d")

# working scripts
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
WORKDIR=${HOMEDIR}/snpEff_impact/${DATASET}/${REF}
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/snpEff_impact/${DATASET}/${REF}/bedfiles
mkdir -p ${WORKDIR}
mkdir -p ${WORKDIR}/logs
mkdir -p ${WORKDIR}/bedfiles

# variables
PASSINGVAR="PASS"
MODIFIERBED=JointCalls_${DATASET}_filterpass_MODIFIER_${CDSTYPE}_snpEff.bed
LOWBED=JointCalls_${DATASET}_filterpass_LOW_${CDSTYPE}_snpEff.bed
MODERATEBED=JointCalls_${DATASET}_filterpass_MODERATE_${CDSTYPE}_snpEff.bed
HIGHBED=JointCalls_${DATASET}_filterpass_HIGH_${CDSTYPE}_snpEff.bed

# output record
LOG=${WORKDIR}/logs/snpEff_impact_tally.txt
echo -e "#IDX\tTotal\tMODIFIER\tLOW\tMODERATE\tHIGH" > ${LOG}

###########################################################
## main
cd ${SCRATCHDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] Merging ${SCRATCHDIR}/JointCalls_${DATASET}_filterpass_impact_${CDSTYPE}_[01-96]_snpEff.bed"

# output and select each types
for ii in $( seq 1 $NJOBS ); do
IDX=$(printf %02d ${ii})
BEDFILE=JointCalls_${DATASET}_filterpass_impact_${CDSTYPE}_${IDX}_snpEff.bed
output_tally
awk 'BEGIN {FS = "\t"}; $4 == "MODIFIER" {print $0}' ${BEDFILE} >> ${MODIFIERBED}
awk 'BEGIN {FS = "\t"}; $4 == "LOW" {print $0}' ${BEDFILE} >> ${LOWBED}
awk 'BEGIN {FS = "\t"}; $4 == "MODERATE" {print $0}' ${BEDFILE} >> ${MODERATEBED}
awk 'BEGIN {FS = "\t"}; $4 == "HIGH" {print $0}' ${BEDFILE} >> ${HIGHBED}
done

# transfer to the WORKDIR
echo -e "[$(date "+%Y-%m-%d %T")] Moving to ${WORKDIR}/bedfiles"
mv -v ${MODIFIERBED} ${WORKDIR}/bedfiles/
mv -v ${LOWBED} ${WORKDIR}/bedfiles/
mv -v ${MODERATEBED} ${WORKDIR}/bedfiles/
mv -v ${HIGHBED} ${WORKDIR}/bedfiles/

cd ${WORKDIR}/bedfiles/

check_length ${MODIFIERBED}
check_length ${LOWBED}
check_length ${MODERATEBED}
check_length ${HIGHBED}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

conda deactivate
