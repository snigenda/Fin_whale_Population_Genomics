# Title: qsub record for process 1
# 
# Author: Meixi Lin
# Date: Tue Nov 19 14:53:25 2019


SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc1/

cd $SCRIPTDIR

bash WGSproc1_qsub_wrapper_20191119.sh GOC038 meixilin Minke
bash WGSproc1_qsub_wrapper_20191119.sh GOC053 meixilin Minke
bash WGSproc1_qsub_wrapper_20191119.sh ENPCA02 meixilin Minke
bash WGSproc1_qsub_wrapper_20191119.sh ENPAK28 meixilin Bryde
bash WGSproc1_qsub_wrapper_20191119.sh GOC038 meixilin Bryde

# Tue Jan  7 15:42:31 2020
bash WGSproc1_qsub_wrapper_20191119.sh ENPCA02 meixilin Bryde

# Wed Jan  8 21:56:33 2020, with new settings 
bash WGSproc1_qsub_wrapper_20200108.sh GOC053 1 meixilin Bryde

# Wed Jan 15 14:45:57 2020
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK28 1 meixilin Minke

# Sun Feb  2 17:53:14 2020
bash WGSproc1_qsub_wrapper_20200108.sh ENPOR10 1 meixilin Minke

# Wed Feb 12 20:04:58 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC002 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK19 1 meixilin Minke

# Thu Feb 13 22:43:19 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC006 1 meixilin Minke

# Fri Feb 21 11:18:43 2020
# Try again for GOC006
bash WGSproc1_qsub_wrapper_20200108.sh GOC006 1 meixilin Minke

# Fri Feb 21 22:09:46 2020
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK20 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK21 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC010 1 meixilin Minke

# Mon Feb 24 10:11:57 2020
# Try update the memory usage for the last step in WGSproc1 
# The original WGSproc1_b was killed at the Picard MergeAlignment step. 
NEXTSCRIPT=WGSproc1_c_MergeAlign_Qualimap_20200216.sh  
qsub -N WGSproc1_c_GOC006 ${NEXTSCRIPT} GOC010 GOC010_A 1 meixilin Minke # Wrong job name! 
qsub -N WGSproc1_c_ENPAK20 ${NEXTSCRIPT} ENPAK20 ENPAK20_A 1 meixilin Minke
qsub -N WGSproc1_c_ENPAK21 ${NEXTSCRIPT} ENPAK21 ENPAK21_A 1 meixilin Minke

# Mon Feb 24 23:38:39 2020
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK22 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC025 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC050 1 meixilin Minke

# Sun Mar  1 15:55:36 2020
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK23 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPAK24 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPCA09 1 meixilin Minke

qalter -l  h_data=24G,h_rt=79200,h_vmem=24G,highp=TRUE,highmem_forced 2227432 

# Wed Mar  4 12:53:43 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC063 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC068 1 meixilin Minke

# Sun Mar  8 15:26:09 2020
# The WGSproc1_c was killed at the Picard MergeAlignment step without error for GOC063. 
cd /u/home/m/meixilin/flashscratch/finwhale/preprocessing/GOC063/Minke
cd temp
rm * 
cd ../
rm 03_c_GOC063_A_MergeBamAlignment.log

SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc1/
cd $SCRIPTDIR

NEXTSCRIPT=WGSproc1_c_MergeAlign_Qualimap_20200216.sh  
qsub -N WGSproc1_c_GOC063 ${NEXTSCRIPT} GOC063 GOC063_A 1 meixilin Minke 


# Sun Mar 29 20:24:43 2020
bash WGSproc1_qsub_wrapper_20200108.sh ENPOR11 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPOR12 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPOR13 1 meixilin Minke

# Wed Apr  1 09:54:28 2020
# the ENPOR11 failed with unknown reasons, resubmit
bash WGSproc1_qsub_wrapper_20200108.sh ENPOR11 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPWA14 1 meixilin Minke

# Thu Apr  2 13:35:32 2020
# ENPOR11 failed becaused of a corrupted gzip file 
# Exception in thread "main" htsjdk.samtools.SAMException: Corrupt GZIP trailer at line 542858621 in fastq /u/project/rwayne/snigenda/finwhale_raw_gnom_seq/FT-SH5556/FT-SA40004-FT-SPN00435_HTJWGDSXX/SH5556_SA40004_S13_L003_R2_001.fastq.gz

bash WGSproc1_qsub_wrapper_20200108.sh ENPWA15 1 meixilin Minke
# Exception in thread "main" htsjdk.samtools.SAMException: invalid code lengths set at line 941516361 in fastq /u/project/rwayne/snigenda/finwhale_raw_gnom_seq/FT-SH5556/FT-SA40008-FT-SPN00435_HTJWGDSXX/SH5556_SA40008_S14_L003_R1_001.fastq.gz

# Thu Apr  2 22:09:48 2020
# check gzip by using the file from a different folder 
bash WGSproc1_qsub_wrapper_20200108.sh GOC112 1 meixilin Minke

# Sun Apr  5 01:26:12 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC116 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC125 1 meixilin Minke

# Fri Apr 10 11:59:40 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC100 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC111 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC091 1 meixilin Minke

# Mon Apr 13 14:47:35 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC082 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC086 1 meixilin Minke

# Tue Apr 14 12:51:45 2020
bash WGSproc1_qsub_wrapper_20200108.sh GOC077 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh GOC080 1 meixilin Minke

# Fri Jun 12 15:58:27 2020
# GOC006 missing the aligned stats (during the hoffman2 change of functions)
# rename the files got transferred to HOMEDIR 
mv GOC006_A_MarkIlluminaAdapters_metrics.txt GOC006_A_MarkIlluminaAdapters_metrics_20200222.txt
mv GOC006_A_Aligned_stats GOC006_A_Aligned_stats_20200222

# resubmit the WGSproc1_a 
NAME="GOC006"
FLAG=0
USER="meixilin"
REF="Minke"
HOMEDIR=/u/project/rwayne/snigenda/finwhale
DICT=${HOMEDIR}/scripts/config/fqpath_fqname_rgid.csv # location of file
WORKSCRIPT=${HOMEDIR}/scripts/WGSproc1/WGSproc1_a_FastqToSam_MarkIlluminaAdapters_20200216.sh
NEXTSCRIPT=${HOMEDIR}/scripts/WGSproc1/WGSproc1_b_Align_20200216.sh 

FQ1=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $3}' $DICT | tr -d \"`
FQ2=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $4}' $DICT | tr -d \"`
RGID=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $5}' $DICT | tr -d \"` 
RGLB=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $6}' $DICT | tr -d \"` 
RGPU=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $7}' $DICT | tr -d \"` 
RGCN=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $8}' $DICT | tr -d \"`
RGPM=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $9}' $DICT | tr -d \"`

qsub -M $USER $WORKSCRIPT ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${USER} ${REF}

# resubmit the WGSproc1_b
# Sun Jun 14 19:11:08 2020
cd /u/scratch/m/meixilin/finwhale/preprocessing/GOC006/Minke
PROGRESSLOG=WGSproc1_${RGID}_${REF}_progress.log

NEXT_JOB_ID=$(qsub -terse -N WGSproc1_b_${NAME} ${NEXTSCRIPT} ${FQ1} ${FQ2} ${NAME} ${RGID} ${RGLB} ${RGPU} ${RGCN} ${RGPM} ${FLAG} ${USER} ${REF}) 
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}

# resubmit the WGSproc1_c
# Tue Jun 16 10:38:50 2020
QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
NAME="GOC006"
RGID=`awk -v pat="$NAME" 'BEGIN {FS = ","}; $2 ~ pat {print $5}' $DICT | tr -d \"` 
FLAG=0
USER="meixilin"
REF="Minke"
NEXTSCRIPT=${HOMEDIR}/scripts/WGSproc1/WGSproc1_c_MergeAlign_Qualimap_20200216.sh
NEXT_JOB_ID=$(${QSUB} -terse -N WGSproc1_c_${NAME} ${NEXTSCRIPT} ${NAME} ${RGID} ${FLAG} ${USER} ${REF}) 
cd /u/scratch/m/meixilin/finwhale/preprocessing/GOC006/Minke
PROGRESSLOG=WGSproc1_${RGID}_${REF}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] Submitted $(basename ${NEXTSCRIPT}): Job ID ${NEXT_JOB_ID}" >> ${PROGRESSLOG}

# transfer the files back 
# Sat Jun 20 23:41:24 2020
cd /u/scratch/m/meixilin/finwhale/preprocessing/GOC006/Minke
mv 01_GOC006_A_FastqToSam.log 01_GOC006_A_FastqToSam_20200616.log
mv 02_GOC006_A_MarkIlluminaAdapters.log 02_GOC006_A_MarkIlluminaAdapters_20200616.log 
mv 03_a_GOC006_A_SamToFastq.log 03_a_GOC006_A_SamToFastq_20200616.log 
mv 03_b_GOC006_A_bwamem.log 03_b_GOC006_A_bwamem_20200616.log
mv 03_c_GOC006_A_MergeBamAlignment.log 03_c_GOC006_A_MergeBamAlignment_20200616.log
mv 04_GOC006_A_qualimap.log 04_GOC006_A_qualimap_20200616.log
mv WGSproc1_GOC006_A_Minke_progress.log WGSproc1_GOC006_A_Minke_progress_20200616.log

cp *.log /u/project/rwayne/snigenda/finwhale/preprocessing/GOC006/Minke/scratchlogs/

###########################################################
# Sat Aug  1 12:23:54 2020
# for the leftover two individuals
SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc1/

cd $SCRIPTDIR
bash WGSproc1_qsub_wrapper_20200108.sh ENPOR11 1 meixilin Minke
bash WGSproc1_qsub_wrapper_20200108.sh ENPWA15 1 meixilin Minke

