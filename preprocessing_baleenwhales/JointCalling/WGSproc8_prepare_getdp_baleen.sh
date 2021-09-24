#! /bin/bash
#$ -l h_rt=23:59:00,h_data=12G,h_vmem=16G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc8_getdp_20210119.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc8_getdp_20210119.err.txt
#$ -m abe

# Usage: generate per contiglist sum on gtDP for all sites
# Please specify the amount of iteration needed based on the DNA_FILE genome
# Example: qsub -t 1-96 WGSproc8_prepare_getdp_baleen.sh

###########################################################
## import packages
sleep $((RANDOM % 120))
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail
###########################################################
## def functions

###########################################################
## def variables
DATASET="f50b4" # 50 fin whales + 4 baleen whales
REF="Minke"
# Variables related to contig list
IDX=$(printf %02d ${SGE_TASK_ID})

# Directories
REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
WORKDIR=${HOMEDIR}/filteredvcf/${DATASET}/${REF}
OUTDIR=${WORKDIR}/variant_summary
LOGDIR=${WORKDIR}/logs
LOG=${LOGDIR}/0x_${DATASET}_${REF}_${BAMHEAD}_getGTDP_${IDX}.log

# Variables
TODAY=$(date "+%Y%m%d")
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)
WORKSCRIPT=${HOMEDIR/baleen_genomes}/scripts/baleen_genomes/JointCalling/extract_vcf_gtdp.py

###########################################################
## main

# logging
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start prepare_getdp before WGSproc8 for ${DATASET} ${REF}; git commit id: ${COMMITID}; Workscript = ${WORKSCRIPT}"

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start prepare_getdp before WGSproc8 for ${DATASET} ${REF}; git commit id: ${COMMITID}; Workscript = ${WORKSCRIPT}" > ${LOG}

mkdir -p ${OUTDIR}
cd ${OUTDIR}

python ${WORKSCRIPT} \
--VCF ${WORKDIR}/JointCalls_${DATASET}_07_snpEff_SIFT_formatted_${IDX}.vcf.gz \
--outfile "${OUTDIR}/${DATASET}_${REF}_chr${IDX}_${TODAY}" \
--contiglist "$IDX" &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi

# also get the mean depth by bcftools stats
bcftools stats -s- ${WORKDIR}/JointCalls_${DATASET}_07_snpEff_SIFT_formatted_${IDX}.vcf.gz --threads 3 > ${DATASET}_${REF}_chr${IDX}_${TODAY}_bcfstats.txt
# grep the per-sample part
grep 'PSC' ${DATASET}_${REF}_chr${IDX}_${TODAY}_bcfstats.txt > ${DATASET}_${REF}_chr${IDX}_${TODAY}_bcfdepth.txt

# cleanup
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done prepare_getdp for ${DATASET} ${REF}"
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done prepare_getdp for ${DATASET} ${REF}" >> ${LOG}

conda deactivate