#!/bin/bash
#$ -l highp,h_data=2G,h_rt=08:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/fsc26/param_btsp/step3_fsc_gather_partialresults_n20rd_20210825.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/fsc26/param_btsp/step3_fsc_gather_partialresults_n20rd_20210825.err.txt
#$ -m abe

# @version      v1
# @script       qsub step3_fsc_gather_partialresults_n20rd_20210825.sh <refid>
# @description  Gather results from replicated runs
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri May 14 00:03:49 2021
# @modification: Fri Jul 16 00:35:29 2021
# @modification: Update for the random start datasets (note also the edits in the R script)
# @modification: Thu Aug 26 07:16:07 2021
# @modification: Not restricting on finishing everything

###########################################################
## import packages

sleep $((RANDOM % 30))

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
REFID=${1}
SETTING='neutral'

# get model settings from the csv file
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
SOURCE0=/u/project/rwayne/pnunez/Results/Demography/fsc2 # the best ll file folder
FSCFILE=${HOMEDIR}/scripts/fsc26/param_btsp/fsc_output_list.csv

if [[ ! -f ${FSCFILE} ]]; then
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: ${FSCFILE} does not exist"
    exit 1
fi

MODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $2}' $FSCFILE)
SUBMODEL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $3}' $FSCFILE)
POP=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $4}' $FSCFILE)
PREFIX=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $5}' $FSCFILE)
DIR=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $6}' $FSCFILE)
BESTLL=$(awk -v pat="$REFID" 'BEGIN {FS = ","}; $1 == pat {print $11}' $FSCFILE)

# define bootstrap prefix
BOOTPREFIX=${SUBMODEL}.${POP}_${PREFIX}

# get original output bestlhoods file
BESTLLFILE=${SOURCE0}/${DIR}/${BESTLL}

if [[ ! -f ${BESTLLFILE} ]]; then
    echo -e "[$(date "+%Y-%m-%d %T")] ERROR: ${BESTLLFILE} does not exist"
    exit 1
fi

# previous simulation settings
NSIM=100 # number of bootstrap SFSes
NREP=20 # number of replicate fsc to check (note the folders are now named as n20 as well)
TODAY=$(date "+%Y%m%d")

SUMDIR=${HOMEDIR}/fsc26/param_btsp/${SETTING}/resultsSummaries/${SUBMODEL}.${POP}/n20rd_${TODAY}
WORKDIR=${HOMEDIR}/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}/${SUBMODEL}.${POP}/${BOOTPREFIX}_bootFSC_n20rd
LOGDIR=${HOMEDIR}/reports/fsc26/param_btsp/step3_fsc_gather_partialresults_n20rd_20210825
mkdir -p ${LOGDIR}

# workscript to gather and plot the results
WORKSCRIPT=${HOMEDIR}/scripts/fsc26/param_btsp/fsc_gather_partialresults_20210823.R
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
LOG=${LOGDIR}/${BOOTPREFIX}_bootFSC_n20rd_Summary_${TODAY}.log

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; REFID = ${REFID}; SETTING = ${SETTING}; WORKDIR = ${WORKDIR}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; REFID = ${REFID}; SETTING = ${SETTING}; WORKDIR = ${WORKDIR}" > ${LOG}


mkdir -p ${SUMDIR}
cd ${SUMDIR}

OUTFILE=${SUMDIR}/${BOOTPREFIX}_bootFSC_n20rd_Summary_${TODAY}.csv

HEADER=$(head -n1 ${WORKDIR}/boot_1/run_1/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot.bestlhoods)
echo -e "bootNum\trunNum\t${HEADER}" | tr "\\t" "," > ${OUTFILE}

# start with the original file
results=$(grep -v [A-Z] ${BESTLLFILE}) # -v invert matching, start lines that does not contains alphabets
echo -e "0\t0\t$results"| tr "\\t" "," >> ${OUTFILE}


for (( ii = 1; ii <= $NSIM; ii++ ))
do
    for (( jj = 1; jj <= $NREP; jj++ ))
    do
        INFILE=${WORKDIR}/boot_${ii}/run_${jj}/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot.bestlhoods
        # only use the results if the input file exists
        if [ -f "$INFILE" ]; then
            results=$(grep -v [A-Z] ${INFILE})
            echo -e "$ii\t$jj\t$results"| tr "\\t" "," >> ${OUTFILE}
        else
            echo -e "ERROR-${ii}-${jj}: ${INFILE} missing"
            echo -e "ERROR-${ii}-${jj}: ${INFILE} missing" >> ${LOG}
            # exit 1 # not exiting because of missing files
        fi
    done
done

wc -l ${OUTFILE} # should be 1 header + 1 estimate + 2000 bootstrap = 2002 lines
wc -l ${OUTFILE} >> ${LOG}

############################################################
# read and estimate the model output using an R script
RINPUT=${BOOTPREFIX}_bootFSC_n20rd_Summary_${TODAY}.csv
Rscript --vanilla ${WORKSCRIPT} ${RINPUT} &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
cat fsc26_CI_*.csv

# ############################################################
# # tarballing the original fsc output directory to reduce file numbers
# # Update: probably not needed for now
# # -c: Create an archive.
# # -z: Compress the archive with gzip.
# # -v: verbose
# # -f: Allows you to specify the filename of the archive.
# cd ${HOMEDIR}/fsc26/param_btsp/${SETTING}/${MODEL}/${POP}/${SUBMODEL}.${POP}/

# tar -czf ${BOOTPREFIX}_bootFSC_n20rd.tar.gz ./${BOOTPREFIX}_bootFSC_n20rd

############################################################
# clean up
conda deactivate
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"



