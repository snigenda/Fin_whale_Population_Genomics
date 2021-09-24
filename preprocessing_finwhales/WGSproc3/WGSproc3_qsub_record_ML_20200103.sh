# Title: qsub record for process 3
# 
# Author: Meixi Lin
# Date: Fri Jan  3 13:18:23 2020

########################################
# Define functions 
########################################

# for resubmitting certain contigs 
# takes argument $1: contig number 
# takes argument $2: individuals 
# default: meixilin Minke
resubmit () {
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc3/WGSproc3_HaplotypeCaller_20200103.sh
cd /u/home/m/meixilin/flashscratch/finwhale/preprocessing/${2}/Minke/
rm 06_${2}_Minke_MarkDuplicates_HaplotypeCaller_${1}.log
rm ${2}_MarkDuplicates_${1}.g.vcf.gz
rm ${2}_MarkDuplicates_${1}.g.vcf.gz.tbi
rm WGSproc3_${2}_Minke_MarkDuplicates_${1}_progress.log
rm WGSproc4_${2}_progress.log
qsub -t ${1} ${WORKSCRIPT} ${2} meixilin Minke
}



SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc3/

cd $SCRIPTDIR

# these were run previously with Other/WGSproc3_HaplotypeCaller_ENPAK28.sh


bash WGSproc3_qsub_wrapper_20200103.sh GOC038 meixilin Minke
bash WGSproc3_qsub_wrapper_20200103.sh GOC038 meixilin Bryde
bash WGSproc3_qsub_wrapper_20200103.sh ENPCA02 meixilin Minke

# Thu Jan  9 10:57:25 2020
bash WGSproc3_qsub_wrapper_20200103.sh ENPCA02 meixilin Bryde

# Fri Jan 10 13:32:48 2020
bash WGSproc3_qsub_wrapper_20200103.sh GOC053 meixilin Minke

# Fri Jan 17 10:51:09 2020
# rerun the contigs that failed 
QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc3/WGSproc3_HaplotypeCaller_20200103.sh

NAME=ENPAK28 ## only need to change this line for submitting different individual's file
USER=meixilin ## "meixilin"/"snigenda"
REF=Minke # should be Minke all the time but for compatibility of the 6 individuals

${QSUB} -t 15 ${WORKSCRIPT} ${NAME} ${USER} ${REF}
${QSUB} -t 41 ${WORKSCRIPT} ${NAME} ${USER} ${REF}

# Thu Feb 20 13:59:39 2020
# Adapt scripts for the hoffman2 memory limits 
# Testing record 
NAME=GOC002
USER=meixilin
REF=Minke
BAMHEAD="MarkDuplicates" # MarkDuplicates/RemoveBadReads

# Conclusion: multiple threading not working very well. and 1/3 of the total resource have to be not reserved. 
# Sat Feb 15 12:50:12 2020
bash WGSproc3_qsub_wrapper_20200103.sh GOC002 meixilin Minke
bash WGSproc3_qsub_wrapper_20200103.sh ENPAK19 meixilin Minke

# for change resources 
qalter -l h_rt=22:00:00,h_data=24G,highp,highmem_forced 2140798
qalter -l h_rt=22:00:00,h_data=24G,highp 2140667

# Wed Feb 26 23:49:13 2020
# GOC010 26
# ENPAK20 65,72 
# Failed without notice 

# GOC010 16,17,21,24,27,34,48,51,74
# GOC006 15,68 
# Failed because of exceeding time limit 

HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc3/WGSproc3_HaplotypeCaller_20200103.sh

qsub -t 65 ${WORKSCRIPT} ENPAK20 meixilin Minke
qsub -t 72 ${WORKSCRIPT} ENPAK20 meixilin Minke

# Thu Feb 27 14:48:03 2020
qsub -t 15 ${WORKSCRIPT} GOC006 meixilin Minke
qsub -t 68 ${WORKSCRIPT} GOC006 meixilin Minke

# check GOC010 jobs remaining time 
# 26 was never finished
test=(16 17 21 24 26 27 34 48 51 74)
cd /u/home/m/meixilin/flashscratch/finwhale/preprocessing/GOC010/Minke/
for ii in ${test[@]}; do
	file=06_GOC010_Minke_MarkDuplicates_HaplotypeCaller_${ii}.log
	tail -1 $file 
done

# output 
# INFO  00:20:53,670 ProgressMeter - NW_006725676.1:3320054      2.0455629E7    22.0 h           64.4 m       93.7%    23.5 h      89.2 m 
# INFO  00:21:16,693 ProgressMeter - NW_006725687.1:23069966              0.0    22.0 h        15250.3 w       90.5%    24.3 h       2.3 h 
# INFO  00:38:30,270 ProgressMeter - NW_006725843.1:2347239       1.806327E7    22.0 h           72.9 m       76.1%    28.9 h       6.9 h 
# INFO  00:42:18,772 ProgressMeter - NW_006725965.1:3511849      2.0890546E7    22.0 h           63.1 m       85.4%    25.7 h       3.8 h 
# INFO  16:42:21,989 ProgressMeter - NW_006726020.1:13521085              0.0    13.9 h        15250.3 w       58.3%    23.9 h      10.0 h 
# INFO  00:46:43,453 ProgressMeter - NW_006726098.1:1725586       2.534341E7    22.0 h           52.0 m       94.2%    23.3 h      81.7 m 
# INFO  01:04:22,200 ProgressMeter - NW_006726353.1:32089884              0.0    22.0 h        15250.3 w      100.0%    22.0 h      15.0 s 
# INFO  01:35:06,148 ProgressMeter - NW_006727209.1:699003      2.8682065E7    22.0 h           46.0 m       98.8%    22.3 h      16.5 m 
# INFO  01:46:29,870 ProgressMeter - NW_006727463.1:35153636              0.0    22.0 h        15250.3 w       97.5%    22.6 h      34.4 m 
# INFO  04:43:05,550 ProgressMeter - NW_006730789.1:43114373              0.0    22.0 h        15250.3 w       88.9%    24.7 h       2.7 h 

# delete the unsuccessful files and resubmit the jobs 
test=(16 17 21 24 26 27 34 48 51 74)
HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc3/WGSproc3_HaplotypeCaller_20200103.sh
cd /u/home/m/meixilin/flashscratch/finwhale/preprocessing/GOC010/Minke/
for ii in ${test[@]}; do
	ls *${ii}*
	rm *${ii}*
	qsub -t ${ii} ${WORKSCRIPT} GOC010 meixilin Minke
done

# Fri Feb 28 11:08:12 2020
# alter the job request to highmem 
jobid=(2206973 2207039 2207040 2207041 2207043 2207044 2207045 2207046 2207047 2207048 2207049 2212418)

for ii in ${jobid[@]}; do
	qalter -l h_data=23G,h_rt=172800,h_vmem=23G,highp=TRUE,highmem_forced $ii
done

# Sun Mar  1 15:31:41 2020
# resubmit the contig 48 for ENPAK22
resubmit 48 ENPAK22 # see function definitions above 

# Wed Apr  1 10:03:00 2020
# resubmit contig 91 for ENPOR12
resubmit 91 ENPOR12 
# new job id: 2650234

# Mon Apr  6 21:54:35 2020
# resubmit the contig 65 for GOC112
resubmit 65 GOC112

# Sat Aug  8 20:36:58 2020
# resubmit 
resubmit 62 ENPOR11
resubmit 82 ENPOR11











