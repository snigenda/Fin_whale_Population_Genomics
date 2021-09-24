#!/bin/bash
#$ -l h_data=8G,h_rt=23:59:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/BUSCO_archive_20210119.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/BUSCO_archive_20210119.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub step3_archive_busco_20210119.sh
# @description	archive the busco analyses into tar balls
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Jan 19 00:19:22 2021

###########################################################
## import packages
set -euo pipefail

###########################################################
## def functions

###########################################################
## def variables
# Don't archive the Minke2 and Blue2 since the output was the same as Minke and Blue
REFLIST=("Minke" "Bryde" "Fin" "Blue" "Humpback")
# BUSCODBLIST=("tetrapoda_odb10" "mammalia_odb10" "laurasiatheria_odb10" "cetartiodactyla_odb10")
BUSCODBLIST=("cetartiodactyla_odb10" "mammalia_odb10")

TODAY=$(date "+%Y%m%d")
REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/busco
WORKDIR=${HOMEDIR}/baleen_genomes/busco

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}"

cd ${SCRATCHDIR}

tar -czf busco_downloads_${TODAY}.tar.gz ./busco_downloads
md5sum busco_downloads_${TODAY}.tar.gz > busco_output_${TODAY}.md5sum

for REF in ${REFLIST[@]}; do
for BUSCODB in ${BUSCODBLIST[@]}; do
tar -czf ${REF}_tran_${BUSCODB:0:5}_${TODAY}.tar.gz ./${REF}_tran_${BUSCODB:0:5}
tar -czf ${REF}_prot_${BUSCODB:0:5}_${TODAY}.tar.gz ./${REF}_prot_${BUSCODB:0:5}
md5sum ${REF}_tran_${BUSCODB:0:5}_${TODAY}.tar.gz >> busco_output_${TODAY}.md5sum
md5sum ${REF}_prot_${BUSCODB:0:5}_${TODAY}.tar.gz >> busco_output_${TODAY}.md5sum
done
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"