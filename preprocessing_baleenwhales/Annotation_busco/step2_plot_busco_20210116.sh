#!/bin/bash
#
# @version 		v0
# @script		step2_plot_archive_busco_20210116.sh
# @description	For plotting and archiving the busco output
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Jan 16 13:33:41 2021
# Ran interactively

###########################################################
## import packages
eval "$(/u/project/rwayne/meixilin/miniconda3/bin/conda shell.bash hook)"
conda activate busco

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
REFLIST=("Minke" "Bryde" "Fin" "Blue" "Humpback" "Minke2" "Blue2")
# BUSCODBLIST=("tetrapoda_odb10" "mammalia_odb10" "laurasiatheria_odb10" "cetartiodactyla_odb10")
BUSCODBLIST=("cetartiodactyla_odb10" "mammalia_odb10")

TODAY=$(date "+%Y%m%d")
REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/busco
WORKDIR=${HOMEDIR}/baleen_genomes/busco/BUSCO_summaries_${TODAY}
mkdir -p ${WORKDIR}

PLOTSCRIPT=/u/project/rwayne/meixilin/miniconda3/envs/busco//bin/generate_plot.py

###########################################################
## main
cd ${WORKDIR}
for REF in ${REFLIST[@]}; do
for BUSCODB in ${BUSCODBLIST[@]}; do
rsync -ahv ${SCRATCHDIR}/${REF}_tran_${BUSCODB:0:5}/short_summary.*_odb10.*.txt ./
rsync -ahv ${SCRATCHDIR}/${REF}_prot_${BUSCODB:0:5}/short_summary.*_odb10.*.txt ./
done
done

python3 ${PLOTSCRIPT} --working_directory ${WORKDIR}
