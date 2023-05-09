#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale/Neutral_regions
#$ -l h_rt=23:00:00,h_data=20G,h_vmem=30G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/NoCpGRepeats_sites.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/NoCpGRepeats_sites.err.txt
#$ -m abe
#$ -t 1-96

# Extract genomic regions that are at least 10 Kb apart from exons and are not within repetitive regions or cpg islands, to new vcf files per chromosome or scaffold
# Adapted from Annabel Beichman's script (for analyzing exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: qsub Extract_noCpG_noRepetitive.sh [reference species]


########## Setting environement

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -oe pipefail


########## Set variables, files and directories

REF=${1}

HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKDIR=${HOMEDIR}/Neutral_regions
VCFDIR=${HOMEDIR}/filteredvcf/all50/${REF}
OUTDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
IDX=$(printf %02d ${SGE_TASK_ID}) 
twentyKb=${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromCDR.0based.bed

mkdir -p ${WORKDIR}
mkdir -p ${OUTDIR}


if [ $REF == 'Minke' ]; then
    REFDIR=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0
    REFERENCE=${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic.fasta
    CPG_REPEATS_ALL=${REFDIR}/CpG_repeats_all.bed
    
fi

if [ $REF == 'Bryde' ]; then
    # Note that for the Bryde's whale, we only used the WMdust.bed file 
    REFDIR=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data
    REFERENCE=${REFDIR}/Balaenoptera_edeni_HiC.fasta
    GFF=${REFDIR}/Balaenoptera_edeni/genome/maker/Balaenoptera_edeni_Balaenoptera_edeni_HiC.fasta_v2.functional.gff3.gz
    MASK=${REFDIR}/Balaenoptera_edeni_HiC_repeats/Balaenoptera_edeni_HiC_repeats_WMdust.bed
    
fi


##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Start extracting nor repetitive regions neither cpg islands for ${REF} ${SGE_TASK_ID} JOB_ID: ${JOB_ID}"
echo "The qsub input"
echo "${REF} ${SGE_TASK_ID}"

PROGRESSLOG=./logs/Extract_No_repeats_cpg_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


########## Obtain vcf files per chromosome or scaffold that do not contain repeat regions or cpg islands and are at least 10 Kb apart from exons

echo -e "[$(date "+%Y-%m-%d %T")]  Extracting neutral regions with GATK SelecVariants... " >> ${PROGRESSLOG}
LOG=./logs/Step1_ExtractNeutralSites_${REF}_SelectVariants_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx15g -Djava.io.tmpdir=./temp -T SelectVariants \
-R ${REFERENCE} \
--variant ${VCFDIR}/JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz \
-XL ${CPG_REPEATS_ALL} \
-L ${twentyKb} \
-o ${OUTDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz &>> ${LOG}


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo "[$(date "+%Y-%m-%d %T")] Done extracting no repeats regions or cpg islands for ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"


conda deactivate



