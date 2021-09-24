SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc7/
# /u/project/rwayne/meixilin/fin_whale/script_not_git

cd $SCRIPTDIR

qsub -t 2 WGSproc7_snpEff_annotations_Minke_20200402.sh Minke meixilin
qsub -t 3-96 WGSproc7_snpEff_annotations_Minke_20200402.sh Minke meixilin

# Sun Apr  5 13:32:42 2020
# rm previous records 
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/testsix/Minke
ls JointCalls_07_snpEff* | wc -l 
rm JointCalls_07_snpEff*
ls JointCalls_0x_A_VariantFiltration* | wc -l 
rm JointCalls_0x_A_VariantFiltration*

cd logs
ls 0x_A_testsix_Minke_MarkDuplicates_variant_filter_* | wc -l 
rm 0x_A_testsix_Minke_MarkDuplicates_variant_filter_*
rm WGSproc7_Minke_MarkDuplicates_*.log

# qsub again 
qsub -t 1-96 WGSproc7_snpEff_annotations_20200402.sh Minke meixilin

# rename Bryde's whale annotations into the same framework 
# check script `mv_bryde_files.sh` 

# the contig 23 were not annotated before 
qsub -t 23 WGSproc7_snpEff_annotations_20200402.sh Bryde meixilin
qalter -l h_data=10G,h_rt=18000,h_vmem=24G 2689758

# Sun Apr  5 20:16:37 2020
# PLEASE USE USER FIRST REF SECOND FROM NOW ON TO MATCH OTHER SCRIPTS 

# Wed May 13 17:49:45 2020 (Done)
qsub -t 1-96 WGSproc7_snpEff_annotations_20200512.sh meixilin Minke

# Wed May 13 17:49:52 2020
# Check job completion status 
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all48/Minke/logs
cat WGSproc7*progress*.log | grep "FAIL"
# The previous jobs completed without a glitch 

###########################################################
# Mon Aug 10 17:55:00 2020
# submit scripts for all50 
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc7/
qsub -t 1-96 WGSproc7_snpEff_annotations_20200810.sh meixilin Minke

# NOTE Some error reported in 
# -rw-r--r-- 1 meixilin rwayne  189 Aug 10 23:25 WGSproc7_all50.err.txt
# Checked, not a problem 

###########################################################
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc7/
qsub -t 96 WGSproc7_snpEff_SIFT_annotations_20200817.sh meixilin Minke
qsub -t 1-96 WGSproc7_snpEff_SIFT_annotations_20200817.sh meixilin Minke

###########################################################
# Tue Aug 18 10:49:12 2020
# the Job ID: 4135479.47 appeared to be glitched
cd /u/scratch/m/meixilin/finwhale/filteredvcf/all50/Minke
rm -r 47/
rm JointCalls_all50_07_snpEff_47.vcf

cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke/logs
rm WGSproc7_Minke_MarkDuplicates_47_progress_all50.log
cd annotation
rm 04_A_all50_Minke_MarkDuplicates_snpEff_47.log  
rm Minke_chr47.genes.txt
rm 04_B_all50_Minke_MarkDuplicates_SIFT_47.log    
rm Minke_chr47.html

cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc7/
qsub -t 47 WGSproc7_snpEff_SIFT_annotations_20200817.sh meixilin Minke

###########################################################
# Wed Aug 19 00:30:22 2020
qsub -t 96 WGSproc7_B_format_SIFT_annotations_20200818.sh meixilin Minke
qsub -t 1-95 WGSproc7_B_format_SIFT_annotations_20200818.sh meixilin Minke

for ii in {01..96}
do 
if ! grep "Done WGSproc7_B for all50 Minke" WGSproc7_Minke_MarkDuplicates_${ii}_progress_all50.log; then
echo "$ii is not done"
fi
done

IDS=("02" "05" "06" "10" "11" "12" "14" "17" "20" "21" "27" "30" "34" "38" "41" "44" "49" "51" "52" "55" "58" "59" "61" "71" "91" "93" "94" "95") 

for ii in ${IDS[@]}
do 
# cat annotation/04_C_all50_Minke_MarkDuplicates_formatSIFT_${ii}.log
# ls JointCalls_all50_07_snpEff_SIFT_formatted_${ii}.vcf.gz
done

for ii in ${IDS[@]}
do 
sleep 2
qsub -t ${ii} WGSproc7_B_format_SIFT_annotations_20200818.sh meixilin Minke
done

qsub -t 20 WGSproc7_B_format_SIFT_annotations_20200818.sh meixilin Minke

###########################################################
###########################################################
# November 2020, rerun files
###########################################################
# Thu Nov  5 13:51:14 2020
# Update the SIFT annotations without the trailing ^M from the output
# First reorganize files
# In hoffman2
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke
ls JointCalls_all50_07_snpEff_SIFT_formatted_*.vcf.gz*
rm JointCalls_all50_07_snpEff_SIFT_formatted_*.vcf.gz*

# move other logs and final output files to archive
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50
mkdir Minke_archive20200824
mv Minke/logs/05_*.log Minke_archive20200824/logs/
mv Minke/logs/0x_*_Minke_get_gtdp.log Minke_archive20200824/logs/
mv Minke/logs/WGSproc6_Minke_MarkDuplicates_*_progress_all50.log Minke_archive20200824/logs/
mv Minke/logs/WGSproc7_Minke_MarkDuplicates_*_progress_all50.log Minke_archive20200824/logs/
mv Minke/logs/WGSproc8_Minke_MarkDuplicates_*_progress_all50.log Minke_archive20200824/logs/
mv Minke/logs/WGSproc9_Minke_MarkDuplicates_archive_sirius_all50.log Minke_archive20200824/logs/
rsync -ahv Minke/variant_summary Minke_archive20200824/
rm -r Minke/variant_summary
rsync -ahv Minke/logs/annotation Minke_archive20200824/logs/
rm -r Minke/logs/annotation
qsub /u/project/rwayne/meixilin/fin_whale/script_not_git/mv_archive_WGSproc8.sh

# Thu Nov  5 14:52:21 2020
cd /u/project/rwayne/snigenda/finwhale/reports
mkdir archive
mv WGSproc9_all50* archive/
mv WGSproc8_all50* archive/
mv WGSproc7_all50* archive/
mv prep_filter_vcf_getdp_all50* archive/

# Thu Nov  5 15:38:24 2020
# check if files transferred back has the same md5sum
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke
grep "JointCalls_all50_06_B_VariantAnnotator_" filteredvcf_all50_Minke_20200824.md5sum | md5sum --check -
mv filteredvcf_all50_Minke_20200824.md5sum ../Minke_archive20200824/

###########################################################
# Fri Nov  6 22:06:23 2020
# testing the script functionalities
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc7
qsub -t 96 WGSproc7_snpEff_SIFT_annotations_format_20201105.sh meixilin Minke

###########################################################
# Fri Nov  6 23:30:42 2020
# final WGSproc7
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc7

qsub -t 1-96 WGSproc7_snpEff_SIFT_annotations_format_20201105.sh meixilin Minke

###########################################################
# Sat Nov  7 20:03:00 2020
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke/logs
for ii in {01..96}
do
if ! grep "Done WGSproc7 for all50 Minke" WGSproc7_Minke_MarkDuplicates_${ii}_progress_all50.log; then
echo "$ii is not done"
fi
done

# 46 was not finished
ii="46"
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke/logs/annotations
rm *46*
cd /u/scratch/m/meixilin/finwhale/filteredvcf/all50/Minke
rm JointCalls_all50_07_snpEff_46.vcf
rm -r 46

cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc7
qsub -t 46 WGSproc7_snpEff_SIFT_annotations_format_20201105.sh meixilin Minke

