#! /bin/bash
#$ -wd /u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH
#$ -l h_rt=90:00:00,h_data=60G,highp,highmem
#$ -N rzoo
#$ -o /u/project/rwayne/pnunez/reports/rzoo_mix10R_3.err.txt
#$ -e /u/project/rwayne/pnunez/reports/rzoo_mix10R_3.out.txt
#$ -m abe
#$ -M pnunez


# This script runs a rscripts that performs runs of homozigosity
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub Step3_roh_rzoo.sh  


#Defining directories ---------------------

workdir=/u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH/

#Main ---------------------------------------

source /u/local/Modules/default/init/modules.sh
module load R/3.6.0

set -o pipefail

#ÂRscript Step3_roh_rzoo_mix10R.R ${workdri}
#Rscript Step3_roh_rzoo_mix15R.R ${workdir}
Rscript Step3_roh_rzoo_mix10R_base3.R ${workdir}
