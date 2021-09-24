#! /bin/bash
#$ -l h_rt=23:59:00,h_data=12G,h_vmem=16G
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc9_20210127.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc9_20210127.err.txt
#$ -m abe

# @version 		v1
# @usage		generate filter statistics
# @description	WGSproc9_f50b4
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Aug 14 10:58:02 2020
# @modification Tue Jan 26 22:05:54 2021
# @modification Update for the baleen_genomes

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
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads
TODAY=$(date "+%Y%m%d")

REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/filteredvcf/${DATASET}/${REF}
WORKDIR=${HOMEDIR}/filteredvcf/${DATASET}/${REF}
STATDIR=${WORKDIR}/filter_stats_${TODAY}
mkdir -p ${WORKDIR}
mkdir -p ${STATDIR}
mkdir -p ${STATDIR}/logs
mkdir -p ${SCRATCHDIR}

# 54 individuals
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125" "BalAcu02" "BalMus01" "EubGla01" "MegNov01")
VCFLIST=${HOMEDIR/baleen_genomes}/scripts/config/baleen_vcffile.txt

IDX=$(printf %02d ${SGE_TASK_ID}) # this is SGE specific array id, essentially the contig list
MYPREFIX="${DATASET}_${REF}_${IDX}_filter_stats"
COMMITID=$(git --git-dir="${HOMEDIR/baleen_genomes}/scripts/.git" --work-tree="${HOMEDIR/baleen_genomes}/scripts" rev-parse master)

PARSESCRIPT="${HOMEDIR/baleen_genomes}/scripts/baleen_genomes/JointCalling/parse_customVCFfilter_summary_20201110.py"
LOG=${STATDIR}/logs/${MYPREFIX}_${TODAY}.log

###########################################################
## main
cd ${SCRATCHDIR}
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc9_a_report_filter_stats for ${DATASET} ${REF}; git commit id: ${COMMITID}; Workscript = ${PARSESCRIPT}"
# echo -e "[$(date "+%Y-%m-%d %T")] VCF names: ${NAMES[@]}"

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Start WGSproc9_a_report_filter_stats for ${DATASET} ${REF}; git commit id: ${COMMITID}; Workscript = ${PARSESCRIPT}" > ${LOG}

# parse and tally the filter statistics
python ${PARSESCRIPT} --dir "${SCRATCHDIR}" --prefix "${MYPREFIX}" --contig "${IDX}" --outdir "${STATDIR}"

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc9_a_report_filter_stats for ${DATASET} ${REF}"
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${IDX}; Done WGSproc9_a_report_filter_stats for ${DATASET} ${REF}" >> ${LOG}

conda deactivate
