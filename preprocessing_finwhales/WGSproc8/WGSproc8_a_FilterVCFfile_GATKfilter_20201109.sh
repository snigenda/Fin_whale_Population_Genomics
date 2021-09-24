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

### Step 8A: Apply mask and hard filters with VariantFiltration

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Start WGSproc8_a_GATKfilter for all50 ${REF}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc8_${REF}_${BAMHEAD}_${IDX}_progress_all50.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${SGE_TASK_ID}; Start WGSproc8_a_GATKfilter for all50 ${REF}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}" >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] GATK Filter variant ... " >> ${PROGRESSLOG}

LOG=./logs/05_A_all50_${REF}_${BAMHEAD}_variant_filter_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

# Mask for:
# 1. CpG sites
# 2. Repeat masked sites in NCBI
# Filter for:
# First set, gatk recommended hard filter.
# https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations
# 1. QD < 2.0: quality by depth (QUAL/DP),
# 2. FS > 60: fisher strand. Phred-scaled probability that there is strand bias at the site. gatk recommended hard filter.
# 3. MQ < 40: root mean square mapping quality over all the reads at the site. Most at 60.
# 4. MQRankSum < -12.5: MappingQualityRankSumTest. A value close to zero is best and indicates little difference between the mapping qualities.
# 5. ReadPosRankSum < -8.0: the u-based z-approximation from the Rank Sum Test for site position within reads. It compares whether the positions of the reference and alternate alleles are different within the reads.
# 6. SOR > 3.0: StrandOddsRatio. Estimate strand bias. Most variants have an SOR value less than 3.

# Second set, gatk recommended hard filter, parameter set for the dataset.
# 7. QUAL < 30.0: quality should be higher than 30

gatk3 -Xmx8g -Djava.io.tmpdir=./temp -T VariantFiltration \
-R ${REFERENCE} \
--logging_level ERROR \
--mask ${MASK} --maskName "FAIL_CpGRep" \
-filter "QD < 2.0" --filterName "FAIL_QD" \
-filter "FS > 60.0" --filterName "FAIL_FS" \
-filter "MQ < 40.0" --filterName "FAIL_MQ" \
-filter "MQRankSum < -12.5" --filterName "FAIL_MQRS" \
-filter "ReadPosRankSum < -8.0" --filterName "FAIL_RPRS" \
-filter "SOR > 3.0" --filterName "FAIL_SOR" \
-filter "QUAL < 30.0" --filterName "FAIL_qual" \
-L ${INTERVAL} \
-V JointCalls_all50_07_snpEff_SIFT_formatted_${IDX}.vcf.gz \
-o JointCalls_all50_08_A_VariantFiltration_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done WGSproc8_a_GATKfilter for all50 ${REF}" >> ${PROGRESSLOG}
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done WGSproc8_a_GATKfilter for all50 ${REF}"

conda deactivate

