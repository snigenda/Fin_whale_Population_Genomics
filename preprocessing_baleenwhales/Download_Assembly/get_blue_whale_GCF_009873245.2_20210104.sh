#!/bin/bash
#$ -l h_data=8G,h_vmem=10G,h_rt=23:00:00
#$ -wd <homedir>
#$ -o <homedir>/reports/get_blue_whale_GCF_009873245.2_20210104.out.txt
#$ -e <homedir>/reports/get_blue_whale_GCF_009873245.2_20210104.err.txt
#$ -m abe

# @version 		v2
# @usage		qsub get_blue_whale_GCF_009873245.2_20210104.sh
# @description	download the blue whale assembly as refseq
# Author: Meixi Lin
# Date: Thu Jul  9 16:47:14 2020
# @modification Sun Sep 13 21:14:22 2020
# @modification modified from download finwhale refgenomes
# @modification: Mon Jan  4 20:15:05 2021
# @modification: Update to download the refseq version

###########################################################
## import packages

conda activate gentools

set -eo pipefail
###########################################################
## def functions

###########################################################
## def variables

QSUB=/u/systems/UGE8.6.4/bin/lx-amd64/qsub
MAKESEQDICT=<homedir>/scripts/baleen_genomes/Preprocessing/make_seqdict_baleen.sh
WORKDIR="<homedir2>/finwhale/cetacean_genomes/blue_whale_genome"
HOMEDIR="<homedir>"

AR_REFERENCE_SEQ=GCF_009873245.2_mBalMus1.pri.v3_genomic.fna.gz  # to archive
REFERENCE_SEQ=GCF_009873245.2_mBalMus1.pri.v3_genomic.fasta

COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"
echo -e "[$(date "+%Y-%m-%d %T")] git commit id: ${COMMITID}"

cd ${WORKDIR}
# download the files
wget -r ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/

# mv the files
mv ./ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/873/245/GCF_009873245.2_mBalMus1.pri.v3/ GCF_009873245.2_mBalMus1.pri.v3/
rm -r ./ftp.ncbi.nlm.nih.gov

# gunzip file
cd GCF_009873245.2_mBalMus1.pri.v3/
cp ${AR_REFERENCE_SEQ} ${AR_REFERENCE_SEQ/.fna.gz/.fasta.gz}
gunzip ${REFERENCE_SEQ}.gz

# generate index for reference sequence submit the next job script
${QSUB} -N make_seqdict_blue_whale ${MAKESEQDICT} "${WORKDIR}/GCF_009873245.2_mBalMus1.pri.v3" "${REFERENCE_SEQ}"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
conda deactivate
