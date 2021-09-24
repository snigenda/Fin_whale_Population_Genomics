#!/bin/bash
#$ -l h_data=4G,h_vmem=4G,h_rt=04:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/gtf_CDS_step1_output_CDS_from_gtf_20201217.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/get_ALLregions_CDS/gtf_CDS_step1_output_CDS_from_gtf_20201217.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub step1_output_CDS_from_gtf_20201217.sh
# @description	generate bedfiles from gtf directly, match up with Vaquita analyses pipeline
# Author: Meixi Lin (meixilin@ucla.edu)
# Thu Dec 17 23:56:52 2020

###########################################################
## import packages

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions
check_length() {
    local FILENAME=${1}
    echo -e "[$(date "+%Y-%m-%d %T")] Getting total length of ${FILENAME}"
    awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${FILENAME}
}

###########################################################
## def variables
DATASET="all50"
REF="Minke"
CDSTYPE="ALLregions"
TODAY=$(date "+%Y%m%d")

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
OUTDIR=${HOMEDIR}/get_ALLregions_CDS//${DATASET}/${REF}/bedfiles
LOGDIR=${OUTDIR}/logs
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

GTFFILE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.gtf

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; output CDS, start_codon and stop_codon from ${GTFFILE}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

awk '{if($3=="CDS" || $3=="start_codon" || $3=="stop_codon") print}' ${GTFFILE} | wc -l

awk 'BEGIN {OFS="\t"}; {if($3=="CDS" || $3=="start_codon" || $3=="stop_codon") print $1,$4-1,$5}' ${GTFFILE} > Minke.CDSstartstop.Coords.FromGTF.0Based.bed

bedtools sort -i Minke.CDSstartstop.Coords.FromGTF.0Based.bed > Minke.CDSstartstop.Coords.FromGTF.0Based.sorted.bed
bedtools merge -i Minke.CDSstartstop.Coords.FromGTF.0Based.sorted.bed > Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.bed

wc -l Minke.CDSstartstop.Coords.FromGTF.0Based.bed

# Check length
check_length Minke.CDSstartstop.Coords.FromGTF.0Based.bed
check_length Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.bed

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
conda deactivate
