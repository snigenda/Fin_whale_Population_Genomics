#! /bin/bash
#$ -l h_rt=23:59:00,h_data=20G,h_vmem=24G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc8_20210122.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc8_20210122.err.txt
#$ -m abe

# Step 8: Variant filtration first on INFO level, then on FORMAT level
# Adapted from Best Practices for GATK3
# Author: Jacquelin Robinson
# Adapted by Meixi & Sergio for use on fin whale project
# Usage: qsub -t 1-<idx> WGSproc8_b_FilterVCFfile_customfilter_baleen.sh
# @modification Tue Nov  3 10:56:46 2020
# @modification 1. Modify to perform GATK variant filter and genotype filters using two scripts 2. Stop backing up to Sirius

sleep $((RANDOM % 120))
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

################################################################################

### Set variables
DATASET="f50b4" # 50 fin whales + 4 baleen whales
REF="Minke"
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/filteredvcf/${DATASET}/${REF}
WORKDIR=${HOMEDIR}/filteredvcf/${DATASET}/${REF}
mkdir -p ${WORKDIR}
mkdir -p ${SCRATCHDIR}

# 54 individuals
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125" "BalAcu02" "BalMus01" "EubGla01" "MegNov01")
VCFLIST=${HOMEDIR/baleen_genomes}/scripts/config/baleen_vcffile.txt

IDX=$(printf %02d ${SGE_TASK_ID})

if [ $REF == 'Minke' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${REFERENCE/GCF_000493695.1_BalAcu1.0_genomic.fasta/contiglist\/BalAcu1.0_genomic.contiglist}_${IDX}.list
    MASK=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/CpG_repeats_all.bed # WM output + RM output + UCSC CpG output
fi

COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)
FILTER_SCRIPT=${HOMEDIR/baleen_genomes}/scripts/baleen_genomes/JointCalling/customVCFfilter_f50b4_20201109.py
MAXDFILE=${HOMEDIR/baleen_genomes}/scripts/baleen_genomes/JointCalling/f50b4_Minke_maxD_20210125.csv

################################################################################

### Step 8B. Apply custom site- and genotype-level filters
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc8_b_FilterVCFfile_customfilter for ${DATASET} ${REF}; git commit id: ${COMMITID}; Workscript = ${FILTER_SCRIPT}; Max depth = ${MAXDFILE}"
echo -e "[$(date "+%Y-%m-%d %T")] VCF names: ${NAMES[@]}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc8_b_customfilter_${DATASET}_${REF}_${BAMHEAD}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${IDX}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] Custom Filter variant ... " >> ${PROGRESSLOG}

LOG=./logs/05_B_${DATASET}_${REF}_${BAMHEAD}_genotype_filter_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

INFILE2=JointCalls_${DATASET}_08_A_VariantFiltration_${IDX}.vcf.gz
OUTFILE2=JointCalls_${DATASET}_08_B_VariantFiltration_${IDX}.vcf.gz
OUTFILTERFILE=${SCRATCHDIR}/${DATASET}_${REF}_${IDX}_filter_stats.tsv

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

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc8_b_customfilter for ${DATASET} ${REF}" >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc8_b_customfilter for ${DATASET} ${REF}"

conda deactivate

