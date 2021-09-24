#!/bin/bash
#$ -l highp,h_data=1G,h_vmem=1G,h_rt=72:00:00
#$ -pe shared 8
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/fsc26/param_btsp/step2_runfsc_btSFS_2DModels_n20rd_20210823.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/fsc26/param_btsp/step2_runfsc_btSFS_2DModels_n20rd_20210823.err.txt
#$ -m abe
#$ -t 1-2000

# @version 		v0
# @usage		qsub step2_runfsc_btSFS_2DModels_n20rd_20210823.sh <model> <submodel> <population> <prefix> <setting>
# @description	Estimate parameters again from bootstrapped SFS generated from the maximum likelihood parameter estimations (n20rd: Each bootstrap estimate 20 times, use random start positions) For use with 2DModels. Split into 400 threads now to get it done faster.
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Aug 23 23:10:54 2021

############################################################
## import packages
sleep $((RANDOM % 120))

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
REFID=${1}
SETTING='neutral'

HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
FSCFILE=${HOMEDIR}/scripts/fsc26/param_btsp/fsc_output_list.csv

if [[ ! -f ${FSCFILE} ]]; then
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: ${FSCFILE} does not exist"
    exit 1
fi

MODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $2}' $FSCFILE)
SUBMODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $3}' $FSCFILE)
POP=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $4}' $FSCFILE)
PREFIX=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $5}' $FSCFILE)
BOOTPREFIX=${SUBMODEL}.${POP}_${PREFIX}

NREP=20
NSIM=1000000
NCORE=8
TODAY=$(date "+%Y%m%d")

# get the SFS to use and the runNumber (customized by how you want to distribute the -t option and the desired number reps)
xx=$((SGE_TASK_ID - 1))
SFSIDX=$((xx % 100 + 1))
ii=$((xx / 100 + 1))

INDIR=${HOMEDIR}/scripts/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}
WORKDIR=${HOMEDIR}/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}/${SUBMODEL}.${POP}/${BOOTPREFIX}_bootFSC_n20rd/boot_${SFSIDX}
mkdir -p ${WORKDIR}
LOGDIR=${HOMEDIR}/reports/fsc26/param_btsp/step2_runfsc_btSFS_2DModels_n20rd/${BOOTPREFIX}
mkdir -p ${LOGDIR}
LOG=${LOGDIR}/${BOOTPREFIX}.${SGE_TASK_ID}.n20rd_${SETTING}_${TODAY}.log

# the bootstrap SFS file
BOOTSFS=${HOMEDIR}/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}/${SUBMODEL}.${POP}/${BOOTPREFIX}_bootSFS/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot_${SFSIDX}/${BOOTPREFIX}_boot_jointMAFpop1_0.obs
if [[ ! -f ${BOOTSFS} ]]; then
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: ${BOOTSFS} does not exist."
    exit 1
fi

COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# the software
FSC26=/u/project/rwayne/software/finwhale/fastsimcoal/fsc26_linux64/fsc26

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}; Using ${BOOTSFS}; Outputting to ${WORKDIR}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; GIT commit id ${COMMITID}; Using ${BOOTSFS}; Outputting to ${WORKDIR}" > ${LOG}


# -n  --numsims 1000      : number of simulations to perform. Also applies for parameter estimation
# -m  --msfs              : computes minor site frequency spectrum
# -M  --maxlhood          : perform parameter estimation by max lhood from SFS
#                           values between iterations
# -L  --numloops 20       : number of loops (ECM cycles) to perform during
#                           lhood maximization. Default is 20
# -c   --cores 1          : number of openMP threads for parameter estimation
#                           (default=1, max=numBatches, use 0 to let openMP choose optimal value)
cd ${WORKDIR}
# if that run exists before, delete that
if [ -d run_${ii} ]; then
    rm -r run_${ii}
fi
mkdir run_${ii}
cd run_${ii}
rsync -a ${BOOTSFS} ./
rsync -a ${INDIR}/${BOOTPREFIX}_boot.est ./
rsync -a ${INDIR}/${BOOTPREFIX}_boot.tpl ./
# rsync -a ${INDIR}/${BOOTPREFIX}_boot.pv ./
${FSC26} -t ${BOOTPREFIX}_boot.tpl -e ${BOOTPREFIX}_boot.est -n ${NSIM} -m -M -L 60 -c ${NCORE} -q &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done. boot_${SFSIDX} run_${ii} ${BOOTPREFIX}" >> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}
