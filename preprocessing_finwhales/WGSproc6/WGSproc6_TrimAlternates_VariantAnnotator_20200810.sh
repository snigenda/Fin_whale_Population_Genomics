#! /bin/bash
#$ -l h_data=15G,h_vmem=16G,h_rt=23:00:00
#$ -wd /u/project/rwayne/snigenda/finwhale
#$ -o /u/project/rwayne/snigenda/finwhale/reports/WGSproc6_all50.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/WGSproc6_all50.err.txt
#$ -m abe

# Step 6: Trim unused alternate alleles and add VariantType and AlleleBalance annotations
# to INFO column
# Adapted from Best Practices for GATK3
# Author: Jacquelin Robinson
# Adapted by Meixi & Sergio for use on fin whale project 

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

# Use the folder in WGSproc5 that store combined and filtered genotype
WORKDIR=${HOMEDIR}/filteredvcf/all50/${REF}
WORKSIRIDIR=${SIRIUSDIR}/filteredvcf/all50/${REF}
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")

mkdir -p ${WORKDIR}
ssh ${USER}@sirius.eeb.ucla.edu "mkdir -p ${WORKSIRIDIR}" 

IDX=$(printf %02d ${SGE_TASK_ID}) # this is SGE specific array id, essentially the contig list 

if [ $REF == 'Minke' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
    INTERVAL=${REFERENCE/GCF_000493695.1_BalAcu1.0_genomic.fasta/contiglist\/BalAcu1.0_genomic.contiglist}_${IDX}.list
fi

if [ $REF == 'Bryde' ]; then
    REFERENCE=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni_HiC.fasta
    INTERVAL=${REFERENCE/Balaenoptera_edeni_HiC.fasta/contiglist\/Bal_edeni_HiC.contiglist}_${IDX}.list
fi

################################################################################

### Trim Alternates 
# echo the input 
echo "[$(date "+%Y-%m-%d %T")] Start WGSproc6 for ${NAMES} ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"

cd ${WORKDIR}
mkdir -p ./logs
mkdir -p ./temp

PROGRESSLOG=./logs/WGSproc6_${REF}_${BAMHEAD}_${IDX}_progress_all50.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] GATK TrimAlternates... " >> ${PROGRESSLOG}

LOG=./logs/02_all50_${REF}_${BAMHEAD}_TrimAlternates_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx8g -Djava.io.tmpdir=./temp -T SelectVariants \
-R ${REFERENCE} \
-trimAlternates \
-L ${INTERVAL} \
-V JointCalls_all50_05_GenotypeGVCFs_${IDX}.vcf.gz \
-o JointCalls_all50_06_A_TrimAlternates_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

################################################################################

### Variant Annotation

# Annotate variant calls with context information
# This tool is designed to annotate variant calls based on their context (as opposed to functional annotation). Various annotation modules are available; see the "Annotation Modules" page linked in the Tool Documentation sidebar for a complete list.

# -A Which annotations to include in variant calls in the output. 
# -G One or more groups of annotations to apply to variant calls. Which groups of annotations to add to the output variant calls.  

echo -e "[$(date "+%Y-%m-%d %T")] GATK VariantAnnotator... " >> ${PROGRESSLOG}
LOG=./logs/03_all50_${REF}_${BAMHEAD}_VariantAnnotator_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx8g -Djava.io.tmpdir=./temp -T VariantAnnotator \
-R ${REFERENCE} \
-G StandardAnnotation \
-A VariantType \
-A AlleleBalance \
-L ${INTERVAL} \
-V JointCalls_all50_06_A_TrimAlternates_${IDX}.vcf.gz \
-o JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

################################################################################

### Move intermediate files to Sirius 

echo -e "[$(date "+%Y-%m-%d %T")] Copying GVCFs files to ${WORKSIRIDIR}... " >> ${PROGRESSLOG}

scp JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz FAIL" >> ${PROGRESSLOG}
    exit 1
fi

scp JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_all50_06_B_VariantAnnotator_${IDX}.vcf.gz.tbi FAIL" >> ${PROGRESSLOG}
    exit 1
fi

rm JointCalls_all50_06_A_TrimAlternates_${IDX}.vcf.gz
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
rm JointCalls_all50_06_A_TrimAlternates_${IDX}.vcf.gz.tbi
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

rm JointCalls_all50_05_GenotypeGVCFs_${IDX}.vcf.gz
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

rm JointCalls_all50_05_GenotypeGVCFs_${IDX}.vcf.gz.tbi
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo "[$(date "+%Y-%m-%d %T")] Done WGSproc6 for ${NAMES} ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"
conda deactivate 

