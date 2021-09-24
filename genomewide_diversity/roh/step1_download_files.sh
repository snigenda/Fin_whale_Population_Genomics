#!/bin/bash
#
# @version 		v0
# @script		bash step1_download_files.sh
# @description	Download ROH files for plotting
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Jan  6 16:22:41 2021

###########################################################
## import packages

###########################################################
## def functions

###########################################################
## def variables

###########################################################
## main
# Thu Jan  7 09:30:01 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/ROH/rohBcftools
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/ROH/rohBcftools/
rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}

LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/ROH/rohZooroh
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH/mix10R_base3/
rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}