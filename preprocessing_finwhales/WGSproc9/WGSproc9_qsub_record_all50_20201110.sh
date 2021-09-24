###########################################################
# Tue Nov 10 19:50:52 2020
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc9
qsub -t 1-96 WGSproc9_a_report_filter_stats_20201110.sh

###########################################################
# Wed Nov 11 10:22:05 2020
# check if previous jobs finished without errors
qacct -j 5234658 | grep 'exit_status' | grep -c 'exit_status  0'
cd /u/project/rwayne/meixilin/fin_whale/analyses/Summary_stats/all50/Minke/filter_stats_20201110/logs/

LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/Summary_stats/all50/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/Summary_stats/all50/Minke/filter_stats_20201110

rsync -ahv --update -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}

###########################################################
# Wed Nov 11 10:34:19 2020
# Performed locally
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts/WGSproc9/WGSproc9_a_plot_sitefilter_stats_all50_20200816.R
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts/WGSproc9/WGSproc9_a_plot_gtfilter_stats_all50_20200816.R







