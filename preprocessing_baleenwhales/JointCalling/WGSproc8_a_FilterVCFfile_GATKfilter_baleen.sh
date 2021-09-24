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
# Usage: qsub -t 1-<idx> WGSproc8_a_FilterVCFfile_GATKfilter_baleen.sh
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

################################################################################

### Step 8A: Apply mask and hard filters with VariantFiltration

# echo the input
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc8_a_FilterVCFfile_GATKfilter for ${DATASET} ${REF}; git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] VCF names: ${NAMES[@]}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc8_a_GATKfilter_${DATASET}_${REF}_${BAMHEAD}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${IDX}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] GATK Filter variant ... " >> ${PROGRESSLOG}

LOG=./logs/05_A_${DATASET}_${REF}_${BAMHEAD}_variant_filter_${IDX}.log
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

gatk3 -Xmx15G -Djava.io.tmpdir=./temp -T VariantFiltration \
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
-V JointCalls_${DATASET}_07_snpEff_SIFT_formatted_${IDX}.vcf.gz \
-o JointCalls_${DATASET}_08_A_VariantFiltration_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc8_a_GATKfilter for ${DATASET} ${REF}" >> ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc8_a_GATKfilter for ${DATASET} ${REF}"

conda deactivate

