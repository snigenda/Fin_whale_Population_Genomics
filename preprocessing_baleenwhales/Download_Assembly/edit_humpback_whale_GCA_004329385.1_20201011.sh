#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/edit_humpback_whale_GCA_004329385.1_20201011.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/edit_humpback_whale_GCA_004329385.1_20201011.err.txt
#$ -m abe

# 
# @version 		v0
# @script		edit_humpback_whale_GCA_004329385.1_20201011.sh	
# @description	edit the directory structure of humpback whales	
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Oct 11 15:24:51 2020

###########################################################
## import packages 
set -eo pipefail
###########################################################
## def functions 

###########################################################
## def variables 

###########################################################
## main 
cd /u/project/rwayne/snigenda/finwhale/cetacean_genomes/humpback_whale_genome
# download the files 
wget -r -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/329/385/GCA_004329385.1_megNov1/ 

# mv the files 
mv ./ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/329/385/GCA_004329385.1_megNov1/ GCA_004329385.1_megNov1/
rm -r ./ftp.ncbi.nlm.nih.gov

# check if the two "fna.gz" files are the same 
md5sum GCA_004329385.1_megNov1_genomic.fna.gz
md5sum GCA_004329385.1_megNov1/GCA_004329385.1_megNov1_genomic.fna.gz

# mv the files (interactive)
rm GCA_004329385.1_megNov1_genomic.fna.gz
ls GCA_004329385.1_megNov1_genomic.*
mv GCA_004329385.1_megNov1_genomic.* GCA_004329385.1_megNov1/
ls megNov1*
mv megNov1* GCA_004329385.1_megNov1/