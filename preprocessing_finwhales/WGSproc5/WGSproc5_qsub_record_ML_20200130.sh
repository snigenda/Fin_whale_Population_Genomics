#!/bin/bash 

# Thu Jan 30 17:00:39 PST 2020
bash WGSproc5_qsub_wrapper_20200130.sh meixilin Minke

# Thu Jan 30 18:28:39 PST 2020
bash WGSproc5_qsub_wrapper_20200130.sh meixilin Bryde

# Fri Jan 31 18:10:09 2020
# The Bryde whale took more than 12 hrs 
HARD_RESOURCE="highp,h_rt=48:00:00,h_data=30G"
QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc5/WGSproc5_GenotypeGVCFs_20200130.sh
USER="meixilin"
REF="Bryde"

${QSUB} -t 17 -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}
${QSUB} -t 20 -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}
${QSUB} -t 21 -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}
${QSUB} -t 22 -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}

# Sun May 10 12:13:59 2020
# MODIFY THE CONFIG FILE WITH THE DATASET INCLUDED 
DICT=./scripts/config/fqpath_fqname_rgid.csv
awk -F, -v OFS=, '{ $10 = "TRUE" }1' ${DICT} > ${DICT/.csv/_1.csv}
# rename the header to all48 
# add FALSE to ENPOR11 ENPWA15

# find names and add to the worker script 
NAMES=`awk -v pat="TRUE" 'BEGIN {FS = ","}; $10 ~ pat {print $2}' $DICT`

# submit the wrapper scripts 
bash WGSproc5_qsub_wrapper_20200510.sh meixilin Minke

# Sun May 10 20:36:26 2020
# Check for completing status 
cd ${HOMEDIR}/filteredvcf/all48/Minke/logs
ls WGSproc5*progress*.log | wc -l 
cat WGSproc5*progress*.log | grep "FAIL" | wc -l

# Some never finished ... 
cat WGSproc5*progress*.log
for ii in {01..96}; do 
	if [ `wc -l WGSproc5_Minke_MarkDuplicates_${ii}_progress_all48.log | awk '{print $1}'` -ne "5" ]; then
		echo WGSproc5_Minke_MarkDuplicates_${ii}_progress_all48.log
	fi
done

for ii in {01..96}; do 
	grep "Done. " 01_all48_Minke_MarkDuplicates_GenotypeGVCFs_${ii}.log
done

# Tue May 12 14:28:32 2020
# resubmit the 75 and 78 jobs 
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all48/Minke/logs
rm *_75* 
rm *_78*
rm JointCalls_05_GenotypeGVCFs_75.vcf.gz*
rm JointCalls_05_GenotypeGVCFs_78.vcf.gz*

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HARD_RESOURCE="h_rt=23:00:00,h_data=23G"
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc5/WGSproc5_GenotypeGVCFs_20200510.sh
USER=meixilin
REF=Minke

${QSUB} -t 75 -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}
${QSUB} -t 78 -l ${HARD_RESOURCE} ${WORKSCRIPT} ${USER} ${REF}

# Tue May 12 15:26:28 2020
# Check for file integrity 
for ii in {01..96}; do 
	echo $ii
	md5sum JointCalls_05_GenotypeGVCFs_${ii}.vcf.gz | ssh meixilin@sirius.eeb.ucla.edu "(cd ${WORKSIRIDIR} && md5sum -c - )"
	md5sum JointCalls_05_GenotypeGVCFs_${ii}.vcf.gz.tbi | ssh meixilin@sirius.eeb.ucla.edu "(cd ${WORKSIRIDIR} && md5sum -c - )"
done

###########################################################
# Sun Aug  9 21:52:26 2020
# submit the jobs 
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc5/

qsub -t 1-96 WGSproc5_GenotypeGVCFs_20200809.sh meixilin Minke

