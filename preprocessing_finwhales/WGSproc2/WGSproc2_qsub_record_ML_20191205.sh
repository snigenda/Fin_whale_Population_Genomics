# Title: qsub record for process 2
# 
# Author: Meixi Lin
# Date: Mon Dec  2 00:03:17 PST 2019


SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc2/

cd $SCRIPTDIR

bash WGSproc2_qsub_wrapper_20191205.sh GOC038 meixilin Minke
bash WGSproc2_qsub_wrapper_20191205.sh GOC053 meixilin Minke
bash WGSproc2_qsub_wrapper_20191205.sh ENPCA02 meixilin Minke
bash WGSproc2_qsub_wrapper_20191205.sh ENPAK28 meixilin Bryde
bash WGSproc2_qsub_wrapper_20191205.sh GOC038 meixilin Bryde

# Wed Jan  8 17:28:27 2020
bash WGSproc2_qsub_wrapper_20200103.sh ENPCA02 meixilin Bryde

# Sat Feb 15 12:50:12 2020
bash WGSproc2_qsub_wrapper_20200103.sh GOC002 1 meixilin Minke
bash WGSproc2_qsub_wrapper_20200103.sh ENPAK19 1 meixilin Minke
