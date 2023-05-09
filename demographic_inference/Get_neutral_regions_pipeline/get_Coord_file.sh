#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale/Neutral_regions
#$ -l h_rt=10:00:00,h_data=4G,h_vmem=10G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/getNeutralCoords.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/getNeutralCoords.err.txt
#$ -m abe
#$ -t 1-96

# Gets neutral coordinates to bed file
# Adapted from Annabel Beichman's script (to analyze exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: qsub get_Coord_file.sh


########## Setting environment

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -oe pipefail

########## Set variables, files and directories

HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKDIR=${HOMEDIR}/Neutral_regions
VCFDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
OUTDIR=${WORKDIR}/neutralVCFs
SCRIPTDIR=${HOMEDIR}/scripts/get_neutral_regions
NeutralCoord_SCRIPT=${SCRIPTDIR}/obtain_noCpG_noRepetitive_coordinates.py
IDX=$(printf %02d ${SGE_TASK_ID})

##### make directories were this information will be stored

mkdir -p ${WORKDIR}/repeatRegions
mkdir -p ${WORKDIR}/get_fasta


##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

echo "[$(date "+%Y-%m-%d %T")] Start creating bed files for ${SGE_TASK_ID} JOB_ID: ${JOB_ID}"
echo "The qsub input"
echo "${SGE_TASK_ID}"

PROGRESSLOG=./logs/create_beds_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


########## Creates a bed file for each vcf file

echo -e "[$(date "+%Y-%m-%d %T")]  creating bed files ... " >> ${PROGRESSLOG}
LOG=./logs/Step2_Creating_bed_files_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

python ${NeutralCoord_SCRIPT} --VCF ${VCFDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz --outfile ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_${IDX}.bed

bedtools merge -d 10 -i ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_${IDX}.bed > ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_mergedMaxDistance10_${IDX}.bed

cat ${WORKDIR}/repeatRegions/min20kb_DistFromCDR_noCpGIsland_noRepeat_mergedMaxDistance10_*.bed | sort -k1,1 -k2,2n > ${WORKDIR}/get_fasta/min20kb_DistFromCDRs_noCpGIsland_noRepeat_mergedMaxDistance10_sorted.bed 


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


conda deactivate

