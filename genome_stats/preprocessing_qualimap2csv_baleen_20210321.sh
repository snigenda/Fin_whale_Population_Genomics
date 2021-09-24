#!/bin/bash
#
# @version 		v1
# @script		preprocessing_qualimap2csv_baleen_20210321.sh
# @description	get qualimap stats into a csv file
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Mar 21 15:03:36 2021

###########################################################
## import packages
set -eo pipefail
###########################################################
## def functions

###########################################################
## def variables
samples=(EubGla01 BalAcu02 BalMus01 MegNov01)

mydate=$(date "+%Y%m%d")
outdir="/u/project/rwayne/meixilin/fin_whale/analyses/Summary_stats/f50b4/Minke/qualimap_${mydate}"
mkdir -p $outdir
outfile=$outdir/qualimap_summary_${mydate}
echo -e "sample\ttotalread\tmappedread\tduprate\tmeancov" > $outfile.raw.txt

###########################################################
## main
cd $outdir
for sample in ${samples[@]}; do
summaryfile=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/preprocessing/${sample}/Minke/${sample}_A_Aligned_stats/genome_results.txt
totalread=`sed -n 's/number of reads = //p' $summaryfile`
mappedread=`sed -n 's/number of mapped reads = //p' $summaryfile`
duprate=`sed -n 's/duplication rate = //p' $summaryfile`
meancov=`sed -n 's/mean coverageData = //p' $summaryfile`
echo -e "$sample\t$totalread\t$mappedread\t$duprate\t$meancov" >> $outfile.raw.txt
done

cat $outfile.raw.txt | awk '{$1=$1;print}' > $outfile.txt

echo -e "[$(date "+%Y-%m-%d %T")] Qualimap summary done"