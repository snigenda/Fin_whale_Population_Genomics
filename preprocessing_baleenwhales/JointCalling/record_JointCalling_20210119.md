# WGSproc5

```bash
# Tue Jan 19 23:01:33 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 baleen_genomes/JointCalling/WGSproc5_GenotypeGVCF_baleen.sh

# Wed Jan 20 11:26:24 2021
qacct -j 6231687 | grep 'exit_status' | grep -c 'exit_status  0'
```

# WGSproc6

```bash
# Wed Jan 20 11:49:28 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 baleen_genomes/JointCalling/WGSproc6_TrimAlternates_VariantAnnotator_baleen.sh

# Thu Jan 21 11:45:24 2021
qacct -j 6235644 | grep 'exit_status' | grep -c 'exit_status  0'
```

# WGSproc7

```bash
# Thu Jan 21 11:46:11 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 baleen_genomes/JointCalling/WGSproc7_snpEff_SIFT_annotations_format_baleen.sh


qacct -j 6243264 | grep 'exit_status' | grep -c 'exit_status  0'

# Fri Jan 22 13:08:32 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 baleen_genomes/JointCalling/WGSproc7_snpEff_SIFT_annotations_format_baleen.sh
qhold 6253020.11-96

# Sat Jan 23 13:20:36 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 19 baleen_genomes/JointCalling/WGSproc7_snpEff_SIFT_annotations_format_baleen.sh
qsub -t 33 baleen_genomes/JointCalling/WGSproc7_snpEff_SIFT_annotations_format_baleen.sh
qsub -t 46 baleen_genomes/JointCalling/WGSproc7_snpEff_SIFT_annotations_format_baleen.sh

# Sun Jan 24 16:40:26 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -l h_data=30G,h_rt=23:59:00,h_vmem=35G -t 46 baleen_genomes/JointCalling/WGSproc7_snpEff_SIFT_annotations_format_baleen.sh
 # 6271546.46

```

# WGSproc8

```bash
# Sun Jan 24 19:01:37 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -t 1-96 baleen_genomes/JointCalling/WGSproc8_prepare_getdp_baleen.sh

# Mon Jan 25 01:28:07 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
OUTDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/filteredvcf/f50b4/Minke/variant_summary

Rscript --vanilla baleen_genomes/JointCalling/parse_vcf_gtdp.R ${OUTDIR}
```

WGSproc8_a

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -t 1-96 baleen_genomes/JointCalling/WGSproc8_a_FilterVCFfile_GATKfilter_baleen.sh

# JOB ID: 6276828
```

WGSproc8_b

**IMPORTANT: The lowest gtdp is bumped up to 29 if it's lower than 29**

* Files changed:

```
"BalAcu02",26
"EubGla01",23
```

```
"BalAcu02",29
"EubGla01",29
```

```bash
# Tue Jan 26 13:20:03 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -t 1-96 baleen_genomes/JointCalling/WGSproc8_b_FilterVCFfile_customfilter_baleen.sh

qacct -j 6291517 | grep 'exit_status' | grep -c 'exit_status  0'
```

# WGSproc9

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -t 1-96 baleen_genomes/JointCalling/WGSproc9_a_report_filter_stats_baleen.sh
qacct -j 6303544 | grep 'exit_status' | grep -c 'exit_status  0'

# move to Summary_stats anyways
# Wed Jan 27 09:34:07 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/Summary_stats/f50b4/Minke/
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/filteredvcf/f50b4/Minke/filter_stats_20210127

rsync -ahv --exclude=logs -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}

# make the plots
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/baleen_genomes/JointCalling/WGSproc9_a_plot_gtfilter_stats_f50b4_20210127.R
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/baleen_genomes/JointCalling/WGSproc9_a_plot_sitefilter_stats_f50b4_20210127.R
```


# Cleanup

```bash
# Tue Jan 26 13:10:24 2021
# Remove JointCalls_05 and JointCalls_06
ls JointCalls_f50b4_05_GenotypeGVCFs_*.vcf.gz* | wc -l
rm JointCalls_f50b4_05_GenotypeGVCFs_*.vcf.gz*

ls JointCalls_f50b4_06_B_VariantAnnotator_*.vcf.gz* | wc -l
rm JointCalls_f50b4_06_B_VariantAnnotator_*.vcf.gz*

ls JointCalls_f50b4_07_snpEff_SIFT_formatted_*.vcf.gz* | wc -l
rm JointCalls_f50b4_07_snpEff_SIFT_formatted_*.vcf.gz*
```

