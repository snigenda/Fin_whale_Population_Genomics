# Fri Jun 26 12:07:00 2020
# Download blue, humpback and the minke used for reference first 
# Megaptera novaeangliae,SRR5665639
# Megaptera novaeangliae,SRR8386009
# Balaenoptera musculus,SRR5665644
# Balaenoptera acutorostrata scammoni,SRR896642

cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -N Blue_whale_SRA baleen_genomes/get_sra_baleen_genomes.sh SRR5665644

###########################################################
# Tue Jul  7 12:45:36 2020
# The jobs failed before and try again with another genome 
qsub -N Humpback_whale_SRA_639 baleen_genomes/get_sra_baleen_genomes.sh SRR5665639
# Use jackie's blue whale sequences 
scp meixilin@sirius.eeb.ucla.edu:/data3/jarobinson/ncbi/public/sra/SRR5665644.sra ${SCRATCHDIR}/
qsub -N Blue_whale_SRA baleen_genomes/SCRATCH/get_sra_baleen_SRR5665644.sh 

###########################################################
# Fri Jul 10 10:47:21 2020
# the jobs kept failing, switching to wget and then fastq-dump
# download the humpback whale 
SRAID="SRR5665639"
SRAPATH="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5665639/SRR5665639.1"
qsub -N Humpback_whale_SRA_639 baleen_genomes/Download_genomes/get_sra_baleen_genomes.sh ${SRAID} ${SRAPATH}

###########################################################
# Fri Jan  8 17:44:52 2021
# start a new way to download sra files
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
SRAID="SRR1802584"
qsub -N ${SRAID} baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh ${SRAID}
qacct -j 6086831

SRAID="SRR5665640"
qsub -N ${SRAID} baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh ${SRAID}
qacct -j 6086874

###########################################################
# Mon Jan 11 09:37:40 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
SRAIDLIST=("SRR926179" "SRR896642" "SRR5665645" "SRR5665646" "SRR935201" "SRR3161874")

for SRAID in ${SRAIDLIST[@]}; do
sleep 120
qsub -N ${SRAID} baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh ${SRAID}
done

###########################################################
# Mon Jan 11 15:46:10 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
SRAIDLIST2=("SRR5665643" "SRR5665641" "SRR5665642" "SRR11097130" "SRR11430498" "SRR8386009")

for SRAID in ${SRAIDLIST2[@]}; do
sleep 120
qsub -N ${SRAID} baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh ${SRAID}
done

# JOBIDS:
cat /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/get_sra_baleen_genomes_20210108* | grep "JOB ID" | sed -e 's/JOB ID \(.*\);/\1/' | cut -d ' ' -f 3 | tr '\n' ' '
JOBIDS=(6086831 6086874 6114410 6114421 6114431 6114440 6114444 6117721 6117825 6117852 6117869 6117918 6119840 6119841)
for jobid in ${JOBIDS[@]}; do
echo ${jobid}
qacct -j $jobid | grep 'exit_status'
done

# SRAIDs
SRAIDLIST=("SRR5665640" "SRR1802584" "SRR926179" "SRR896642" "SRR5665645" "SRR5665646" "SRR935201" "SRR3161874" "SRR5665643" "SRR5665641" "SRR5665642" "SRR11097130" "SRR11430498" "SRR8386009")

for SRAID in ${SRAIDLIST[@]}; do
cat /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/get_sra_baleen_genomes_20210108* | grep -c ${SRAID}
echo $SRAID
done

###########################################################
# Mon Jan 11 20:37:56 2021
# some were not included
qsub -N SRR3161874 baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh SRR3161874 # JOBID: 6119840
qsub -N SRR8386009 baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh SRR8386009 # JOBID: 6119841

###########################################################
# Tue Jan 12 09:14:51 2021
# one job did not complete: JOBID: 6117869; SRA ID: SRR11097130
# SRR1802584
qsub -N SRR1802584 baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh SRR1802584 # JOBID: 6143034

# Tue Jan 12 15:14:58 2021
# downloaded the wrong one
# SRR11097130
qsub -N SRR11097130 baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh SRR11097130 # JOBID: 6143034

###########################################################
# Tue Jan 12 16:05:57 2021
# download the other two already downloaded from ENA as well
SRAID="SRR5665639"
qsub -N ${SRAID} baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh ${SRAID} # 6156420

SRAID="SRR5665644"
qsub -N ${SRAID} baleen_genomes/Download_SRA/get_sra_baleen_genomes.sh ${SRAID} # 6156421

###########################################################
# Tue Jan 12 23:11:41 2021
# check that all files have been downloaded
SRAIDLIST=("SRR926179" "SRR1802584" "SRR896642" "SRR5665645" "SRR5665646" "SRR5665644" "SRR935201" "SRR3161874" "SRR5665643" "SRR5665641" "SRR5665642" "SRR5665640" "SRR11097130" "SRR11430498" "SRR5665639" "SRR8386009")

for SRAID in ${SRAIDLIST[@]}; do
ls ${SRAID}* | wc -l
done
