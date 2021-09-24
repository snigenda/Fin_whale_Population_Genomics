#! /bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc5_all50.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc5_all50.err.txt
#$ -m abe

# Step 5: GenotypeGVCFs
# Adapted from Best Practices for GATK3
# Generates a multi-sample joint VCF file with genotypes at ALL sites
# Author: Jacquelin Robinson
# Adapted by Meixi & Sergio for use on fin whale project 
# Usage: ./WGSproc5_GenotypeGVCFs.sh [chromosome]

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

# QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub

################################################################################

### Set variables
USER=${1}
REF=${2}
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

HOMEDIR=/u/project/rwayne/snigenda/finwhale
SIRIUSDIR=/data3/finwhale 

# Start a new folder that store combined and filtered genotype
WORKDIR=${HOMEDIR}/filteredvcf/all50/${REF}
WORKSIRIDIR=${SIRIUSDIR}/filteredvcf/all50/${REF}
# DICT=${HOMEDIR}/scripts/config/fqpath_fqname_rgid.csv
# NAMES=`awk -v pat="TRUE" 'BEGIN {FS = ","}; $10 ~ pat {print $2}' $DICT`
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")

mkdir -p ${WORKDIR}
ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${WORKSIRIDIR}" 

IDX=$(printf %02d ${SGE_TASK_ID})

if [ $REF == 'Minke' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${REFERENCE/GCF_000493695.1_BalAcu1.0_genomic.fasta/contiglist\/BalAcu1.0_genomic.contiglist}_${IDX}.list
fi

if [ $REF == 'Bryde' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
    INTERVAL=${REFERENCE/Balaenoptera_edeni_HiC.fasta/contiglist\/Bal_edeni_HiC.contiglist}_${IDX}.list
fi

################################################################################

### Genotype GVCF
# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc5 for ${NAMES} ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"
echo "The qsub input: ${USER} ${REF} ${JOB_ID}.${SGE_TASK_ID}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc5_${REF}_${BAMHEAD}_${IDX}_progress_all50.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] GATK GenotypeGVCFs... " >> ${PROGRESSLOG}

# generate a list of vcf file pathes 

VCFFILES=()
for ii in ${NAMES[@]}; do
	VCFFILES+=(${HOMEDIR}/preprocessing/${ii}/${REF}/${ii}_${BAMHEAD}_${IDX}.g.vcf.gz)
done 

LOG=./logs/01_all50_${REF}_${BAMHEAD}_GenotypeGVCFs_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

# sleep random times to avoid overload gatk3? 
sleep $((RANDOM % 120))
# start genotype GVCF
gatk3 -Xmx8g -Djava.io.tmpdir=./temp -T GenotypeGVCFs \
-R ${REFERENCE} \
-allSites \
-stand_call_conf 0 \
-L ${INTERVAL} \
$(for ff in ${VCFFILES[@]}; do echo "-V ${ff} "; done) \
-o JointCalls_all50_05_GenotypeGVCFs_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

### Save jointVCF files
echo -e "[$(date "+%Y-%m-%d %T")] Copying GenotypeGVCFs files to ${WORKSIRIDIR}... " >> ${PROGRESSLOG}

scp JointCalls_all50_05_GenotypeGVCFs_${IDX}.vcf.gz ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

scp JointCalls_all50_05_GenotypeGVCFs_${IDX}.vcf.gz.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


echo "[$(date "+%Y-%m-%d %T")] Done WGSproc5 for ${NAMES} ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"
conda deactivate 


