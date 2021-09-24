###########################################################
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
###########################################################
# Fri Jul 10 14:09:12 2020
# download the fin whale reference genome 
qsub testsix_comparisons/Baphy_NCBI/step1_download_seq_makedict.sh

###########################################################
# Sun Sep 13 22:41:02 2020
qsub baleen_genomes/Download_Assembly/get_blue_whale_GCA_009873245.2_20200913.sh

###########################################################
# Wed Sep 16 12:12:30 2020
Rscript --vanilla baleen_genomes/Download_Assembly/contiglist_blue_whale_GCA_009873245.2_20200915.R
Rscript --vanilla baleen_genomes/Download_Assembly/contiglist_humpback_whale_GCA_009873245.2_20200915.R # 73

###########################################################
# Mon Jan  4 20:21:57 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub baleen_genomes/Download_Assembly/get_blue_whale_GCF_009873245.2_20210104.sh

###########################################################
# Tue Jan 12 23:38:38 2021
Rscript --vanilla baleen_genomes/Download_Assembly/contiglist_blue_whale_GCF_009873245.2_20210113.R
echo $?
