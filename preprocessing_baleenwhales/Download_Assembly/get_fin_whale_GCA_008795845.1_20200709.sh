#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/testsix/download_finwhale_refgenome_20200709.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/testsix/download_finwhale_refgenome_20200709.err.txt
#$ -m abe

# @version 		v0
# @usage		
# @description	testsix/download_finwhale_refgenome
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Jul  9 16:47:14 2020

###########################################################
## import packages 
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail
###########################################################
## def functions 

###########################################################
## def variables 
WORKDIR="/u/project/rwayne/snigenda/finwhale/cetacean_genomes/fin_whale_genome"
mkdir -p ${WORKDIR}

AR_REFERENCE_SEQ=GCA_008795845.1_Baphy_genomic.fna.gz  # to archive 
REFERENCE_SEQ=GCA_008795845.1_Baphy_genomic.fasta

###########################################################
## main 
cd ${WORKDIR}
# download the files 
wget -r -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/795/845/GCA_008795845.1_Baphy/

# mv the files 
mv ./ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/795/845/GCA_008795845.1_Baphy/ GCA_008795845.1_Baphy/
rm -r ./ftp.ncbi.nlm.nih.gov

# gunzip file
cd GCA_008795845.1_Baphy/
cp ${AR_REFERENCE_SEQ} ${AR_REFERENCE_SEQ/.fna.gz/.fasta.gz}
gunzip ${AR_REFERENCE_SEQ/.fna.gz/.fasta.gz}

# generate index for reference sequence
bwa index ${REFERENCE_SEQ}

samtools faidx ${REFERENCE_SEQ} 

# generate dictionary (step 4c)
picard -Xmx4G CreateSequenceDictionary \
REFERENCE=${REFERENCE_SEQ} \
OUTPUT=${REFERENCE_SEQ/fasta/dict}

conda deactivate