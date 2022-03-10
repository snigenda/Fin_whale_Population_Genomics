# The final files sources (in the raw data)

```bash
/u/project/rwayne/meixilin/fin_whale/analyses/slim/finwhale_2pop_ancestralChange_040121
/u/project/rwayne/meixilin/fin_whale/analyses/slim/finwhale_2pop_ancestralChange_noMig_040121
/u/project/rwayne/meixilin/fin_whale/analyses/slim/finwhale_ENP_recovery_OLDchr_071421
/u/project/rwayne/meixilin/fin_whale/analyses/slim/finwhale_ENP_OLDchr_071421
```

# step1: summarize the output

```R
source('~/Lab/finwhale_manuscript/scripts/slim/summary_output_20210910.R', echo=TRUE)
```

# step2: plot the output

```R
# Figure 5
source('/Users/linmeixi/Lab/finwhale_manuscript/scripts/maintext_plots/Fig5_slim_enp20gen_ancsplit_20210910.R', echo = TRUE)
# Figure S17
source('/Users/linmeixi/Lab/finwhale_manuscript/scripts/slim/FigS17_slim_ancsplit_allelecount_20210912.R', echo = TRUE)
```

# step3: generate the summary statistics for main text

```bash
# Sun Sep 12 14:17:55 PDT 2021
LOG='./data/slim/derived_data/logs/step3_generate_summary_stats_maintext_20210912.log'

Rscript --vanilla /Users/linmeixi/Lab/finwhale_manuscript/scripts/slim/step3_generate_summary_stats_maintext_20210912.R &> ${LOG}
echo $?

```




