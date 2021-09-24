#!/bin/bash 

# Fri Jan 31 18:47:52 2020
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc6/
bash WGSproc6_qsub_wrapper_20200131.sh meixilin Minke

# Sun Feb  2 17:29:27 2020
# Previous code did not copy the vcf.gz.tbi files 
REF="Minke"
HOMEDIR=/u/project/rwayne/snigenda/finwhale
SIRIUSDIR=/data3/finwhale 
WORKDIR=${HOMEDIR}/filteredvcf/testsix/${REF}
WORKSIRIDIR=${SIRIUSDIR}/filteredvcf/testsix/${REF}

cd $WORKDIR
for ii in {1..96}; do
	IDX=$(printf %02d ${ii})
	# scp JointCalls_05_GenotypeGVCFs_${IDX}.vcf.gz.tbi meixilin@sirius.eeb.ucla.edu:${WORKSIRIDIR}
	scp JointCalls_06_B_VariantAnnotator_${IDX}.vcf.gz.tbi meixilin@sirius.eeb.ucla.edu:${WORKSIRIDIR}
	rm JointCalls_05_GenotypeGVCFs_${IDX}.vcf.gz.tbi
	rm JointCalls_06_A_TrimAlternates_${IDX}.vcf.gz.tbi
done

REF="Bryde"
HOMEDIR=/u/project/rwayne/snigenda/finwhale
SIRIUSDIR=/data3/finwhale 
WORKDIR=${HOMEDIR}/filteredvcf/testsix/${REF}
WORKSIRIDIR=${SIRIUSDIR}/filteredvcf/testsix/${REF}

cd $WORKDIR
for ii in {1..23}; do
	IDX=$(printf %02d ${ii})
	scp JointCalls_05_GenotypeGVCFs_${IDX}.vcf.gz.tbi meixilin@sirius.eeb.ucla.edu:${WORKSIRIDIR}
done

# Sun Feb  2 17:46:51 2020
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc6/
bash WGSproc6_qsub_wrapper_20200131.sh meixilin Bryde

# Wed May 13 00:24:49 2020
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc6/
bash WGSproc6_qsub_wrapper_20200510.sh meixilin Minke

# Wed May 13 09:53:05 2020
# check for completion status 
cd ${HOMEDIR}/filteredvcf/all48/Minke/logs
ls WGSproc6*progress*.log | wc -l 
cat WGSproc6*progress*.log | grep "FAIL"

###########################################################
# Mon Aug 10 12:07:31 2020
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc6/

qsub -t 1-96 WGSproc6_TrimAlternates_VariantAnnotator_20200810.sh meixilin Minke
