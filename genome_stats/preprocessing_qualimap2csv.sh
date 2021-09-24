#!/bin/bash
#
# @version 		v1
# @script		preprocessing_qualimap2csv.sh
# @description	get qualimap stats into a csv file
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jun 12 15:37:40 2020
# @modification Fri Aug 28 13:23:01 2020
# @modification Add two individuals

###########################################################
## import packages
set -eo pipefail
###########################################################
## def functions

###########################################################
## def variables
samples=(ENPAK19 ENPAK20 ENPAK21 ENPAK22 ENPAK23 ENPAK24 ENPAK25 ENPAK26 ENPAK27 ENPAK28 ENPAK29 ENPAK30 ENPBC16 ENPBC17 ENPBC18 ENPCA01 ENPCA02 ENPCA03 ENPCA04 ENPCA05 ENPCA06 ENPCA07 ENPCA08 ENPCA09 ENPOR10 ENPOR11 ENPOR12 ENPOR13 ENPWA14 ENPWA15 GOC002 GOC006 GOC010 GOC025 GOC038 GOC050 GOC053 GOC063 GOC068 GOC071 GOC077 GOC080 GOC082 GOC086 GOC091 GOC100 GOC111 GOC112 GOC116 GOC125)
mydate=$(date "+%Y%m%d")
outdir="/u/project/rwayne/meixilin/fin_whale/analyses/Summary_stats/all50/Minke/qualimap_${mydate}"
mkdir -p $outdir
outfile=$outdir/qualimap_summary_${mydate}
echo -e "sample\ttotalread\tmappedread\tduprate\tmeancov" > $outfile.raw.txt

###########################################################
## main
cd $outdir
for sample in ${samples[@]}; do
summaryfile=/u/project/rwayne/snigenda/finwhale/preprocessing/${sample}/Minke/${sample}_A_Aligned_stats/genome_results.txt
totalread=`sed -n 's/number of reads = //p' $summaryfile`
mappedread=`sed -n 's/number of mapped reads = //p' $summaryfile`
duprate=`sed -n 's/duplication rate = //p' $summaryfile`
meancov=`sed -n 's/mean coverageData = //p' $summaryfile`
echo -e "$sample\t$totalread\t$mappedread\t$duprate\t$meancov" >> $outfile.raw.txt
done

cat $outfile.raw.txt | awk '{$1=$1;print}' > $outfile.txt

echo -e "[$(date "+%Y-%m-%d %T")] Qualimap summary done"