# step0: combine the bcftools output to a more refined text output

```bash
# Wed Jan  6 17:32:41 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub important_results/Runs_of_homozygosity/step0_combine_bcftools_20210104.sh

qacct -j 6073011
```

# step1: download the remote data for local processing

# step2: format and categorize the ROH objects for better plotting

```R
source('~/Lab/finwhale_manuscript/scripts/roh/step2_format_categorize_roh_20210105.R', echo=TRUE)
```
# step3: summarize the step2 output, perform wilcoxon ranksum test

```R
source('~/Lab/finwhale_manuscript/scripts/roh/step3_summarize_wilcoxtest_froh_20210325.R', echo=TRUE)
```

# step4: plot ROH concordances in two softwares

```bash
TODAY=$(date "+%Y%m%d")
LOG=/Users/linmeixi/google_drive/finwhale/analyses/important_results/Runs_of_homozygosity/logs/step4_roh_bcfzoo_compare_${TODAY}.log

Rscript --vanilla "/Users/linmeixi/Lab/fin_whale/scripts_analyses/important_results/Runs_of_homozygosity/step4_roh_bcfzoo_compare_20210325.R" &> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
```

# step5: get individual's Froh tables

Using the three categories and two softwares

```R
# Fri Aug 13 16:55:44 2021
source('~/Lab/fin_whale/scripts_analyses/important_results/Runs_of_homozygosity/step7_individual_roh_20210813.R', echo=TRUE)

```

# step6: Figure S8 plot ROH genomewide distribution 

```R
# Sun Sep 12 15:39:46 2021
source('~/Lab/finwhale_manuscript/scripts/roh/step6_FigS8_plot_roh_distribution_20210325.R', echo=TRUE)

```
