#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale/Neutral_regions
#$ -l h_rt=23:00:00,h_data=20G,h_vmem=30G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/GetNeutral_sites.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/GetNeutral_sites.err.txt
#$ -t 1-96
#$ -m abe

# Extract the neutral regions that are not conserved (i.e. do not map to zebra fish genome), to new vcf files per chromosome or scaffold.
# The extracted neutral regions will be used to buil the SFS projection. 
# Adapted from Annabel Beichman's script (to analyze exomes) by Sergio Nigenda to analyze whole genome data.
# Usage in hoffman2: qsub neutral_Sites2vcf.sh [reference species]

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -oe pipefail

########## Set variables, files and directories

REF=${1}

HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKDIR=${HOMEDIR}/Neutral_regions
VCFDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
OUTDIR=${WORKDIR}/neutralVCFs
SCRIPTDIR=${HOMEDIR}/scripts/get_neutral_regions
NeutralCoord_SCRIPT=${SCRIPTDIR}/obtain_noCpG_noRepetitive_coordinates.py
twentyKb=${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromCDR.0based.bed
IDX=$(printf %02d ${SGE_TASK_ID})

if [ $REF == 'Minke' ]; then
    REFDIR=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0
    REFERENCE=${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic.fasta
    
fi

if [ $REF == 'Bryde' ]; then
    # Note that for the Bryde's whale, we only used the WMdust.bed file 
    REFDIR=${HOMEDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data
    REFERENCE=${REFDIR}/Balaenoptera_edeni_HiC.fasta
    GFF=${REFDIR}/Balaenoptera_edeni/genome/maker/Balaenoptera_edeni_Balaenoptera_edeni_HiC.fasta_v2.functional.gff3.gz
    MASK=${REFDIR}/Balaenoptera_edeni_HiC_repeats/Balaenoptera_edeni_HiC_repeats_WMdust.bed
    
fi


##### make directories were this information will be stored

mkdir -p ${WORKDIR}/repeatRegions
mkdir -p ${WORKDIR}/get_fasta
mkdir -p ${WORKDIR}/zebra_fish
mkdir -p ${WORKDIR}/HQ_neutral_sites # This directory will have neutral regions going through 3 checks: CpG Islands, repetitive regions, and do not blast to fish
mkdir -p ${OUTDIR}/beds

Fishy=${WORKDIR}/zebra_fish/fish.matches.eval.1e-10.0based.bed


##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Start extracting no conserved regions for ${REF} ${SGE_TASK_ID} JOB_ID: ${JOB_ID}"
echo "The qsub input"
echo "${REF} ${SGE_TASK_ID}"

PROGRESSLOG=./logs/Extract_No_conserved_regions_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


########## Extracts regions that do not blasted to zebra fish genome

echo -e "[$(date "+%Y-%m-%d %T")]  Extracting final neutral regions with GATK SelecVariants... " >> ${PROGRESSLOG}
LOG=./logs/Step4_ExtractNeutralSites_${REF}_SelectVariants_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx15g -Djava.io.tmpdir=./temp -T SelectVariants \
-R ${REFERENCE} \
--variant ${VCFDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz  \
-XL ${Fishy} \
-L ${twentyKb} \
-o ${OUTDIR}/Neutral_sites_SFS_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


########## Optional step. Creates a bed file containing all the Neutral sites and calculates the total length of the neutral sequences obtained with this pipeline.

python ${NeutralCoord_SCRIPT} --VCF ${OUTDIR}/Neutral_sites_SFS_${IDX}.vcf.gz --outfile ${OUTDIR}/beds/Neutral_sites_SFS_${IDX}.bed

bedtools merge -i ${OUTDIR}/beds/Neutral_sites_SFS_${IDX}.bed > ${OUTDIR}/beds/Neutral_sites_SFS_merge_${IDX}.bed

cat ${OUTDIR}/beds/Neutral_sites_SFS_merge_*.bed | sort -k1,1 -k2,2n > ${OUTDIR}/beds/Neutral_sites_SFS_sorted.bed

cp ${OUTDIR}/beds/Neutral_sites_SFS_sorted.bed ${WORKDIR}/HQ_neutral_sites

cd ${WORKDIR}/HQ_neutral_sites

mv Neutral_sites_SFS_sorted.bed Final_HQ_neutral_regions.bed

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${WORKDIR}/HQ_neutral_sites/Final_HQ_neutral_regions.bed > ${WORKDIR}/HQ_neutral_sites/totalPassingSequence.txt


conda deactivate

