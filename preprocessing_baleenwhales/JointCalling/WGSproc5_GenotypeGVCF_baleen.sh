#! /bin/bash
#$ -l h_rt=23:59:00,h_data=20G,h_vmem=24G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc5_20210119.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc5_20210119.err.txt
#$ -m abe

# Step 5: GenotypeGVCFs for baleen genomes
# Adapted from Best Practices for GATK3
# Generates a multi-sample joint VCF file with genotypes at ALL sites
# Author: Jacquelin Robinson
# Adapted by Meixi & Sergio for use on fin whale project
# Usage: qsub -t 1-<idx> WGSproc5_GenotypeGVCF_baleen.sh

sleep $((RANDOM % 120))
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

# QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub

################################################################################

### Set variables
DATASET="f50b4" # 50 fin whales + 4 baleen whales
REF="Minke"
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes
WORKDIR=${HOMEDIR}/filteredvcf/${DATASET}/${REF}
mkdir -p ${WORKDIR}

# 54 individuals
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125" "BalAcu02" "BalMus01" "EubGla01" "MegNov01")
VCFLIST=${HOMEDIR/baleen_genomes}/scripts/config/baleen_vcffile.txt

IDX=$(printf %02d ${SGE_TASK_ID})

if [ $REF == 'Minke' ]; then
    REFERENCE=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${REFERENCE/GCF_000493695.1_BalAcu1.0_genomic.fasta/contiglist\/BalAcu1.0_genomic.contiglist}_${IDX}.list
fi

COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

################################################################################

### Genotype GVCF
# echo the input
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc5_GenotypeGVCF for ${DATASET} ${REF}; git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] VCF names: ${NAMES[@]}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc5_GenotypeGVCF_${DATASET}_${REF}_${BAMHEAD}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${IDX}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] GATK GenotypeGVCFs... " >> ${PROGRESSLOG}

# generate a list of vcf file pathes
VCFFILES=()
while read ii; do
    VCFFILES+=(${ii}_${IDX}.g.vcf.gz)
done < ${VCFLIST}

LOG=./logs/01_${DATASET}_${REF}_${BAMHEAD}_GenotypeGVCFs_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

# start genotype GVCF
gatk3 -Xmx15G -Djava.io.tmpdir=./temp -T GenotypeGVCFs \
-R ${REFERENCE} \
-allSites \
-stand_call_conf 0 \
-L ${INTERVAL} \
$(for ff in ${VCFFILES[@]}; do echo "-V ${ff} "; done) \
-o JointCalls_${DATASET}_05_GenotypeGVCFs_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

###########################################################
# ### Save jointVCF files
# echo -e "[$(date "+%Y-%m-%d %T")] Copying GenotypeGVCFs files to ${WORKSIRIDIR}... " >> ${PROGRESSLOG}

# scp JointCalls_${DATASET}_05_GenotypeGVCFs_${IDX}.vcf.gz ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# scp JointCalls_${DATASET}_05_GenotypeGVCFs_${IDX}.vcf.gz.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}

# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi
# echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}
###########################################################

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc5_GenotypeGVCF for ${DATASET} ${REF}"
conda deactivate


