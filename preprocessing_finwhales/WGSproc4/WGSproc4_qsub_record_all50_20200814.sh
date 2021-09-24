#!/bin/bash
# 
# @version 		
# @script		WGSproc4_qsub_record_all50_20200814.sh
# @description	record the final archive log for preprocessing files 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Aug 14 14:44:39 2020

###########################################################
## def variables 

SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc4/
cd $SCRIPTDIR

###########################################################
## main
NAMES=("ENPAK19" "ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")

for NAME in ${NAMES[@]}
do
sleep 2
qsub -N WGSproc4_all50_${NAME} WGSproc4_Sirius_Archive_20200814.sh ${NAME} meixilin Minke
done

###########################################################
# Fri Aug 14 16:59:45 2020
# Fix more permission problems 
NAMES=("ENPAK25" "ENPAK26" "ENPAK27" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "GOC071")
echo ${#NAMES[@]} # get variable length 
for NAME in ${NAMES[@]}
do
cd /u/project/rwayne/snigenda/finwhale/preprocessing/${NAME}/Minke
cat ./scratchlogs/WGSproc4_b_${NAME}_Sirius_Archive_progress.log
cat ./scratchlogs/07_${NAME}_${REF}_Sirius_Archive.log
done
# ENPAK19 ENPAK25 GOC038 successed, all other failed 

# delete the line that failed 
NAMES=("ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")
for NAME in ${NAMES[@]}
do
cd /u/project/rwayne/snigenda/finwhale/preprocessing/${NAME}/Minke
sed -i '/FAIL/d' ./scratchlogs/WGSproc4_b_${NAME}_Sirius_Archive_progress.log
done


for NAME in ${NAMES[@]}
do
echo $NAME
cd /u/project/rwayne/snigenda/finwhale/preprocessing/${NAME}/Minke
sed -i '/Broken pipe/d' ./scratchlogs/07_${NAME}_${REF}_Sirius_Archive.log
tail -5 ./scratchlogs/07_${NAME}_${REF}_Sirius_Archive.log
done

###########################################################
# Sat Aug 15 00:54:30 2020
# in sirius 
NAMES=("ENPAK20" "ENPAK21" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK26" "ENPAK27" "ENPAK28" "ENPAK29" "ENPAK30" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA02" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA06" "ENPCA07" "ENPCA08" "ENPCA09" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC002" "GOC006" "GOC010" "GOC025" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")

for NAME in ${NAMES[@]}
do
sleep 1
screen -dmS preprocessing_${NAME} bash -c "bash /home/meixilin/finwhale/scripts/WGSproc4/WGSproc4_Sirius_Archive_check_20200814.sh ${NAME}"
done

###########################################################
# Sat Aug 15 19:13:27 2020
# for brydes 
qsub -N WGSproc4_all50_ENPCA08_Bryde WGSproc4_Sirius_Archive_20200814.sh ENPCA08 meixilin Bryde
qsub -N WGSproc4_all50_GOC071_Bryde WGSproc4_Sirius_Archive_20200814.sh GOC071 meixilin Bryde

cd /u/project/rwayne/snigenda/finwhale/preprocessing/ENPCA08/Bryde

# in sirius check the file move process 
NAMES=("ENPCA08" "GOC071")
for NAME in ${NAMES[@]}
do
sleep 1
screen -dmS preprocessing_${NAME} bash -c "bash /home/meixilin/scripts/WGSproc4_Sirius_Archive_check_20200814.sh ${NAME}"
done
###########################################################
# Sun Aug 16 12:05:31 2020
# delete the remove bad reads folders 
cd /u/project/rwayne/snigenda/finwhale/preprocessing/GOC053/Minke
rm ./GOC053_RemoveBadReads.bam
rm ./GOC053_RemoveBadReads.bai

###########################################################
# Thu Aug 20 00:31:42 2020
# delete some bam files not needed 
NAMES=("ENPAK25" "ENPAK26" "ENPAK27" "ENPAK29" "ENPCA03" "ENPCA04" "ENPCA05" "GOC025" "GOC050" "GOC063")

for NAME in ${NAMES[@]}
do
cd /u/project/rwayne/snigenda/finwhale/preprocessing/${NAME}/Minke
rm ${NAME}_MarkDuplicates.bam
rm ${NAME}_MarkDuplicates.bai
done

###########################################################

NAMES=("ENPAK25" "ENPAK26" "ENPAK27" "ENPAK29" "ENPCA03" "ENPCA04" "ENPCA05" "GOC025" "GOC050" "GOC063")

for NAME in ${NAMES[@]}
do
cd /u/project/rwayne/snigenda/finwhale/preprocessing/${NAME}/Minke
rm ${NAME}_MarkDuplicates.bam
rm ${NAME}_MarkDuplicates.bai
done

###########################################################
# Fri Sep 18 00:16:58 2020
# delete all the bam files not used in following analyses
cd /u/project/rwayne/snigenda/finwhale/preprocessing
NAMES=("ENPAK19" "ENPAK20" "ENPAK22" "ENPAK23" "ENPAK24" "ENPAK25" "ENPAK26" "ENPAK27" "ENPAK29" "ENPBC16" "ENPBC17" "ENPBC18" "ENPCA01" "ENPCA03" "ENPCA04" "ENPCA05" "ENPCA07" "ENPOR10" "ENPOR11" "ENPOR12" "ENPOR13" "ENPWA14" "ENPWA15" "GOC006" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC111" "GOC116")
echo ${#NAMES[@]}

for NAME in ${NAMES[@]}
do
cd /u/project/rwayne/snigenda/finwhale/preprocessing/${NAME}/Minke
rm ${NAME}_MarkDuplicates.bam
rm ${NAME}_MarkDuplicates.bai
done


cd /u/project/rwayne/snigenda/finwhale/preprocessing
rm -r ENPCA08/Bryde/
rm -r GOC071/Bryde/

###########################################################
