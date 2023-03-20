#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=23:00:00
#$ -wd <homedir>
#$ -o <homedir>/reports/get_blue_whale_GCA_009873245.2_20200913.out.txt
#$ -e <homedir>/reports/get_blue_whale_GCA_009873245.2_20200913.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub get_blue_whale_GCA_009873245.2_20200913.sh
# @description	download the blue whale assembly with annotations
# Author: Meixi Lin
# Date: Thu Jul  9 16:47:14 2020
# @modification Sun Sep 13 21:14:22 2020
# @modification modified from download finwhale refgenomes

###########################################################
## import packages

conda activate gentools

set -eo pipefail
###########################################################
## def functions

###########################################################
## def variables
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
MAKESEQDICT=<homedir>/scripts/baleen_genomes/Preprocessing/make_seqdict_baleen.sh
WORKDIR="<homedir2>/finwhale/cetacean_genomes/blue_whale_genome"

AR_REFERENCE_SEQ=GCA_009873245.2_mBalMus1.v2_genomic.fna.gz  # to archive
REFERENCE_SEQ=GCA_009873245.2_mBalMus1.v2_genomic.fasta

###########################################################
## main
cd ${WORKDIR}
# download the files
wget -r -nv ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/873/245/GCA_009873245.2_mBalMus1.v2/

# mv the files
mv ./ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/873/245/GCA_009873245.2_mBalMus1.v2/ GCA_009873245.2_mBalMus1.v2/
rm -r ./ftp.ncbi.nlm.nih.gov

# gunzip file
cd GCA_009873245.2_mBalMus1.v2/
cp ${AR_REFERENCE_SEQ} ${AR_REFERENCE_SEQ/.fna.gz/.fasta.gz}
gunzip ${REFERENCE_SEQ}.gz

# generate index for reference sequence submit the next job script
${QSUB} -N make_seqdict_blue_whale ${MAKESEQDICT} "${WORKDIR}/GCA_009873245.2_mBalMus1.v2" "${REFERENCE_SEQ}"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
conda deactivate
