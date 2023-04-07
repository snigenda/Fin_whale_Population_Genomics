#! /bin/bash

#$ -wd /u/project/rwayne/snigenda/finwhale/Neutral_regions
#$ -l h_data=10G,h_vmem=16G,h_rt=23:00:00
#$ -o /u/project/rwayne/snigenda/finwhale/reports/GetHiQcoords.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/GetHiQcoords.err.txt
#$ -t 1-96
#$ -m abe

# Author: Sergio Nigenda, based on Tanya Phung's scripts for step 10 of her NGS pipeline 
# Generates a bed file with the high quality sites in the VCF files
# Usage example: qsub GetHiQualCoords_20200602.sh Minke


########## Setting conda environment 

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

# QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub


########## Set variables, directories and files

REF=${1}

HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKDIR=${HOMEDIR}/Neutral_regions/HiQualCoords
SCRIPTDIR=${HOMEDIR}/scripts/get_neutral_regions
VCFDIR=${HOMEDIR}/filteredvcf/all50/${REF}
HiQualCoords_SCRIPT=${SCRIPTDIR}/obtain_high_qual_coordinates.py
IDX=$(printf %02d ${SGE_TASK_ID})

mkdir -p ${WORKDIR}


########## Logging

# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start GetHiQCoords for ${REF} {IDX} Job ID: ${SGE_TASK_ID}"
echo "The qsub input"
echo "${REF} ${SGE_TASK_ID}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/GetHiQualSitesCoords_A_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${SGE_TASK_ID}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] Selecting high quality coordinates ... " >> ${PROGRESSLOG}

LOG=./logs/01_A_Get_HighQuality_Coords_${REF}_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}


########## Step 9A. Selecting high quality coordinates and missing sites. If ${HiQualCoords_SCRIPT} is defined as obtain_high_qual_coordinates.py it selects only high quality data to later build neutral SFS without projection. However, if is defined as obtain_high_qual_coordinates_miss.py then it includes high quality and missing data for SFS projection. 

python ${HiQualCoords_SCRIPT} --VCF ${VCFDIR}/JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz --outfile HQsitesCoords_${IDX}.bed


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


########## Step 9B. Merging High quality coordinates (The previous step is made for each variant, therefore here we merge those coordinates to have less individual regions)

PROGRESSLOG=./logs/GetHiQualSitesCoords_B_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] Merging high quality coordinates ... " >> ${PROGRESSLOG}

LOG=./logs/01_B_Merge_HighQuality_Coords_${REF}_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

bedtools merge -i HQsitesCoords_${IDX}.bed > HQsitesCoords_merged_${IDX}.bed

wait

cat HQsitesCoords_merged_*.bed | sort -k1,1 -k2,2n > all_HQCoords_sorted_merged.bed

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

conda deactivate
