#! /bin/bash
#$ -l h_rt=23:59:00,h_data=20G,h_vmem=24G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc7_20210122.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc7_20210122.err.txt
#$ -m abe

# Step 7: Add annotations for the mutations using snpEff and SIFT
# Adapted from Best Practices for GATK3
# Authors: Sergio & Meixi
# Usage: qsub -t 1-<idx> WGSproc7_snpEff_SIFT_annotations_format_baleen.sh

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

# software directories
SNPEFFDIR=/u/project/rwayne/software/finwhale/miniconda2/envs/gentools/share/snpeff-4.3.1t-2
SIFT=/u/project/rwayne/software/finwhale/sift/SIFT4G_Annotator.jar

# define variables by references
if [ $REF == 'Minke' ]; then
    SNPEFFDB=Baac01.10776
    SIFTDB=/u/project/rwayne/jarobins/utils/programs/sift/databases/GCF_000493695.1_BalAcu1.0
fi

# workscript
WORKSCRIPT=${HOMEDIR/baleen_genomes}/scripts/baleen_genomes/JointCalling/formatSIFTvcf_f50b4_20201105.py

# commit id
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

################################################################################

### logging
# echo the input
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc7_snpEff_SIFT_annotations_format for ${DATASET} ${REF}; git commit id: ${COMMITID}"
echo -e "[$(date "+%Y-%m-%d %T")] VCF names: ${NAMES[@]}"

cd ${WORKDIR}

# starts logging
LOGDIR=${WORKDIR}/logs/annotations
mkdir -p ${LOGDIR}
PROGRESSLOG=${WORKDIR}/logs/WGSproc7_snpEff_SIFT_${DATASET}_${REF}_${BAMHEAD}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}.${IDX}" > ${PROGRESSLOG}
echo -e "[$(date "+%Y-%m-%d %T")] snpEff annotations based on ${SNPEFFDB} ... " >> ${PROGRESSLOG}

LOG=${LOGDIR}/04_A_${DATASET}_${REF}_${BAMHEAD}_snpEff_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

################################################################################

### snpEFF annotations

snpEff -Xmx15G -nodownload -v -canon -stats ${LOGDIR}/${REF}_chr${IDX}.html \
${SNPEFFDB} \
JointCalls_${DATASET}_06_B_VariantAnnotator_${IDX}.vcf.gz > \
${SCRATCHDIR}/JointCalls_${DATASET}_07_snpEff_${IDX}.vcf 2> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

################################################################################

### SIFT annotations
LOG=${LOGDIR}/04_B_${DATASET}_${REF}_${BAMHEAD}_SIFT_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] SIFT annotations based on ${SIFTDB} ... " >> ${PROGRESSLOG}

mkdir -p ${SCRATCHDIR}/${IDX}
cd ${SCRATCHDIR}/${IDX}

java -Xmx15G -jar ${SIFT} \
-c -t -r ./ -d ${SIFTDB} -i ${SCRATCHDIR}/JointCalls_${DATASET}_07_snpEff_${IDX}.vcf &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi
mv JointCalls_${DATASET}_07_snpEff_${IDX}_SIFTpredictions.vcf JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.vcf
rm ${SCRATCHDIR}/JointCalls_${DATASET}_07_snpEff_${IDX}.vcf

################################################################################

### format SIFT annotations
echo -e "[$(date "+%Y-%m-%d %T")] Formatting JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.vcf using ${WORKSCRIPT} ... " >> ${PROGRESSLOG}

LOG=${LOGDIR}/04_C_${DATASET}_${REF}_${BAMHEAD}_formatSIFT_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

python ${WORKSCRIPT} JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.vcf | \
bgzip > JointCalls_${DATASET}_07_snpEff_SIFT_formatted_${IDX}.vcf.gz
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

tabix -p vcf JointCalls_${DATASET}_07_snpEff_SIFT_formatted_${IDX}.vcf.gz
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}
date "+%Y-%m-%d %T" >> ${LOG}

# delete files
echo -e "[$(date "+%Y-%m-%d %T")] Cleaning up ... " >> ${PROGRESSLOG}
set -eo pipefail
rm JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.vcf

# move files and logs
mv Events.log ${LOGDIR}/Events_SIFTannotations_${IDX}.log
mv JointCalls_${DATASET}_07_snpEff_${IDX}_SIFTannotations.xls ${LOGDIR}/JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.xls
mv JointCalls_${DATASET}_07_snpEff_SIFT_formatted_${IDX}.vcf.gz* ${WORKDIR}/

echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}
date "+%Y-%m-%d %T" >> ${LOG}

# ################################################################################

# ### Move intermediate files to Sirius
# cd ${WORKDIR}

# echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_${DATASET}_07_snpEff_SIFT_*.vcf.gz files to ${WORKSIRIDIR}... " >> ${PROGRESSLOG}

# scp JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.vcf.gz ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_${DATASET}_07_snpEff_${IDX}.vcf.gz FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# scp JointCalls_${DATASET}_07_snpEff_SIFT_${IDX}.vcf.gz.tbi ${USER}@sirius.eeb.ucla.edu:${WORKSIRIDIR}
# exitVal=${?}
# if [ ${exitVal} -ne 0 ]; then
#     echo -e "[$(date "+%Y-%m-%d %T")] Copying JointCalls_${DATASET}_07_snpEff_${IDX}.vcf.gz.tbi FAIL" >> ${PROGRESSLOG}
#     exit 1
# fi

# # rm JointCalls_${DATASET}_06_B_VariantAnnotator_${IDX}.vcf.gz
# # rm JointCalls_${DATASET}_06_B_VariantAnnotator_${IDX}.vcf.gz.tbi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc7_snpEff_SIFT_annotations_format for ${DATASET} ${REF}"
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc7_snpEff_SIFT_annotations_format for ${DATASET} ${REF}" >> ${PROGRESSLOG}

conda deactivate
