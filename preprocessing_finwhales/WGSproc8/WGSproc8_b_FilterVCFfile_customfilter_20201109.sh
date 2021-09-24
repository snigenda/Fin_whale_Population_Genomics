#! /bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc8_all50.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc8_all50.err.txt
#$ -m abe

# Step 8: Variant filtration first on INFO level, then on FORMAT level
# Adapted from Best Practices for GATK3
# Author: Jacquelin Robinson
# Adapted by Meixi & Sergio for use on fin whale project
# Usage: ./WGSproc8_FilterVCF.sh [chromosome]
# @modification Tue Nov  3 10:56:46 2020
# @modification 1. Modify to perform GATK variant filter and genotype filters using two scripts 2. Stop backing up to Sirius

sleep $((RANDOM % 120))
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail


################################################################################

### Set variables
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")

USER=${1}
REF=${2}
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

HOMEDIR=/u/project/rwayne/snigenda/finwhale
# SIRIUSDIR=/data3/finwhale

# Use the folder in WGSproc5 that store combined and filtered genotype
WORKDIR=${HOMEDIR}/filteredvcf/all50/${REF}
# WORKSIRIDIR=${SIRIUSDIR}/filteredvcf/all50/${REF}
SCRATCHDIR=/u/scratch/${USER:0:1}/${USER}/finwhale/filteredvcf/all50/${REF}

mkdir -p ${WORKDIR}
mkdir -p ${SCRATCHDIR}
# ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${WORKSIRIDIR}"

IDX=$(printf %02d ${SGE_TASK_ID}) # this is SGE specific array id, essentially the contig list
COMMITID=`git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master`

if [ $REF == 'Minke' ]; then
    MINKEDIR=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0
    REFERENCE=${MINKEDIR}/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${MINKEDIR}/contiglist/BalAcu1.0_genomic.contiglist_${IDX}.list
    MASK=${MINKEDIR}/CpG_repeats_all.bed # WM output + RM output + UCSC CpG output
fi

if [ $REF == 'Bryde' ]; then
    # Note that for the Bryde's whale, we only used the WMdust.bed file
    BRYDEDIR=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data
    REFERENCE=${BRYDEDIR}/Balaenoptera_edeni_HiC.fasta
    INTERVAL=${BRYDEDIR}/contiglist/Bal_edeni_HiC.contiglist_${IDX}.list
    MASK=${BRYDEDIR}/Balaenoptera_edeni_HiC_repeats/Balaenoptera_edeni_HiC_repeats_WMdust.bed
fi

################################################################################

### Step 8B. Apply custom site- and genotype-level filters
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Start WGSproc8_b_customfilter for all50 ${REF}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc8_${REF}_${BAMHEAD}_${IDX}_progress_all50.log
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Start WGSproc8_b_customfilter for all50 ${REF}" >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] Custom Filter variant ... " >> ${PROGRESSLOG}

LOG=./logs/05_B_all50_${REF}_${BAMHEAD}_genotype_filter_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

FILTER_SCRIPT=${HOMEDIR}/scripts/WGSproc8/customVCFfilter_all50_20201109.py

INFILE2=JointCalls_all50_08_A_VariantFiltration_${IDX}.vcf.gz
OUTFILE2=JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz
MAXDFILE=${HOMEDIR}/scripts/WGSproc8/maxD_file/maxD_custom_filter_${REF}_all50.csv
OUTFILTERFILE=${SCRATCHDIR}/all50_${REF}_${IDX}_filter_stats.tsv

# perform custom filter
python ${FILTER_SCRIPT} ${INFILE2} ${MAXDFILE} ${OUTFILTERFILE} | bgzip > ${OUTFILE2}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] Custom Filter FAIL" >> ${PROGRESSLOG}
    exit 1
fi

tabix -p vcf ${OUTFILE2}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] tabix FAIL" >> ${PROGRESSLOG}
    exit 1
fi
date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo -e "[$(date "+%Y-%m-%d %T")] Get summary statistics ... " >> ${PROGRESSLOG}
# check the filter stats and the custom filter output are the same
echo ${OUTFILTERFILE} >> ${LOG}
tail -1 ${OUTFILTERFILE} >> ${LOG}
wc -l ${OUTFILTERFILE} >> ${LOG}
# gzip the output filter stats files
gzip ${OUTFILTERFILE}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] gzip FAIL" >> ${PROGRESSLOG}
    exit 1
fi

# get simple summary of the final output vcf
echo ${OUTFILE2} >> ${LOG}
bcftools +counts ${OUTFILE2} &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] bcftools counts FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

# ################################################################################

# ### Move intermediate files to Sirius

# echo -e "[$(date "+%Y-%m-%d %T")] Copying filtered files to ${WORKSIRIDIR}... " >> ${PROGRESSLOG}

# scp ${INFILE2} ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying ${INFILE2} FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi
# scp ${INFILE2}.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying ${INFILE2}.tbi FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi
# scp ${OUTFILE2} ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying ${OUTFILE2} FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# scp ${OUTFILE2}.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying ${OUTFILE2}.tbi FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# # rm JointCalls_all50_07_snpEff_${IDX}.vcf.gz
# # rm JointCalls_all50_07_snpEff_${IDX}.vcf.gz.tbi
# # rm ${INFILE2}
# # rm ${INFILE2}.tbi

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done WGSproc8_b_customfilter for all50 ${REF}" >> ${PROGRESSLOG}
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done WGSproc8_b_customfilter for all50 ${REF}"

conda deactivate

