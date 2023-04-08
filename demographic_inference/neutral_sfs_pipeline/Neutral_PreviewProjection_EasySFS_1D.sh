#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale/SFS
#$ -l h_rt=20:00:00,h_data=10G,h_vmem=15G
#$ -N PreviewProjection
#$ -o /u/project/rwayne/snigenda/finwhale/reports/Preview_Neutral_Projection.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/Preview_Neutral_Projection.err.txt
#$ -t 1-96
#$ -m abe

# Usage: qsub easySFS_1_ProjectionPreview_onlyPASSpopHet75.sh
# This script will get and parse project preview for a given vcf file 
# @modification Mon Jul 13 21:44:56 2020
# @modification Final settings. Add per population maxHet filters at 0.75.

set -eo pipefail

###########################################################
## import packages 
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

###########################################################
## define functions
 
# select variants 
# $1 = input VCF name
# $2 = output VCF name
gatk_select_variants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, select SNPs for ${1}" >> ${PROGRESSLOG}
gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude_sample_name ${ss} "; done) \
--excludeFiltered \
--removeUnusedAlternates \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
--variant ${1} -o ${2} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Finish gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, select SNPs for ${1}" >> ${PROGRESSLOG}
}

# select invariant sites
# $1 = input VCF name
# $2 = output VCF name
gatk_select_invariants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, Select Invariant sites for ${1}" >> ${PROGRESSLOG}
gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude_sample_name ${ss} "; done) \
--excludeFiltered \
--removeUnusedAlternates \
--selectTypeToInclude NO_VARIATION \
--selectTypeToExclude INDEL \
--variant ${1} -o ${2} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Finishh gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, Select Invariant sites for ${1}" >> ${PROGRESSLOG}
}

# perform the projection; parse the projection and plot the projection
# $1 = input VCF name
# $2 = output file suffix (OUTDIRSNPS previously defined)
easySFS_projection_popHet75() {
# generate easySFS preview 
# -a keep all snps 
# -v verbose 
echo -e "[$(date "+%Y-%m-%d %T")] Start easySFS, perPop maxHet 0.75 for ${1}" >> ${PROGRESSLOG}
python ${WORKSCRIPT} -i ${1} -maxHetFilter ${maxHetFilter} -p ${POPFILE} --preview -a -v > ${OUTDIRSNPS}/PreviewProjection_${2}.txt 2>> ${LOG}

# parse the preview (the population had to match presupplied variables)
INFILE=${OUTDIRSNPS}/PreviewProjection_${2}.txt
echo -e "[$(date "+%Y-%m-%d %T")] Parsing easySFS output ${INFILE}" >> ${PROGRESSLOG}
for pop in ${POPS[@]}
do
OUTFILE=${OUTDIRSNPS}/${pop}_${2}_PreviewProjection_Rformat.txt
echo "projection,snps" > ${OUTFILE}
grep -A1 "$pop$" ${INFILE} | tail -1 | \
sed 's/(//g' | sed 's/)//g' | sed 's/, /,/g' |  tr '\t' '\n' >> ${OUTFILE}
done

echo -e "[$(date "+%Y-%m-%d %T")] Finish easySFS, perPop maxHet 0.75 for ${1}" >> ${PROGRESSLOG}

# perform plotting (redirect stderr to log as well)
# echo -e "[$(date "+%Y-%m-%d %T")] plotting ${INFILE}" >> ${PROGRESSLOG}
# cd ${OUTDIRSNPS}
# Rscript --vanilla ${PLOTSCRIPT} "${OUTDIRSNPS}" "${OUTDIRSNPS}/OptimalProjectionPlots/" ${2} >> ${LOG} 2>&1
}

###########################################################
## def variables 
HOMEDIR=/u/project/rwayne/snigenda/finwhale
VCFDIR=${HOMEDIR}/Neutral_regions/neutralVCFs
WORKDIR=${HOMEDIR}/SFS
WORKSCRIPT=${HOMEDIR}/scripts/SFS/easySFS_a.py
PLOTSCRIPT=${HOMEDIR}/scripts/SFS/easySFS_function_determineOptimalProjection.R
IDX=$(printf %02d ${SGE_TASK_ID})

# Minke whale reference genome 
REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
EXCLUDE_SAMPLE=("ENPCA01" "ENPCA09" "ENPOR12" "GOC010" "GOC080" "GOC111")
POPS=("ENP" "GOC")

OUTDIRSNPS=${WORKDIR}/Neutral/SNPs
OUTDIRINVS=${WORKDIR}/Neutral/Invariants
mkdir -p ${OUTDIRSNPS}
mkdir -p ${OUTDIRINVS}
mkdir -p ${WORKDIR}/Neutral/logs

POPFILE=${HOMEDIR}/scripts/config/pop_map.txt
maxHetFilter=0.75 # still need max het filter if using only PASS sites in filtered vcf, because of population specific het (popHet) could still exceed maxHet 

# additional arguments 
mydate=$(date "+%Y%m%d")
LOG=${WORKDIR}/Neutral/logs/easySFS_ProjectionPreview_onlyPASSpopHet75_${mydate}_${IDX}.log
PROGRESSLOG=${WORKDIR}/Neutral/logs/ProjectionPreview_${mydate}_${IDX}_progress.log
###########################################################
## main 
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; easySFS projection Final (onlyPASS, popHet=0.75)" 

cd ${OUTDIRSNPS}
gatk_select_variants_onlyPASS ${VCFDIR}/Neutral_sites_SFS_${IDX}.vcf.gz ${OUTDIRSNPS}/SNPs_neutral_for_SFS_${IDX}.vcf.gz 

easySFS_projection_popHet75 ${OUTDIRSNPS}/SNPs_neutral_for_SFS_${IDX}.vcf.gz "neutral_SNPS_${IDX}"

cd ${OUTDIRINVS}
gatk_select_invariants_onlyPASS ${VCFDIR}/Neutral_sites_SFS_${IDX}.vcf.gz ${OUTDIRINVS}/Invariants_neutral_for_SFS_${IDX}.vcf.gz

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} ${IDX} Done" 
echo -e "[$(date "+%Y-%m-%d %T")] job for ${JOB_ID} ${IDX} Done" >> ${PROGRESSLOG}

conda deactivate 
