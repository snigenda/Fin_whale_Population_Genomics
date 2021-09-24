Author: Meixi Lin
Last updated: Tue Aug 24 14:50:10 2021

# Important settings
1. Window size: 1Mb
2. Step size: 1Mb
3. Maximum missing: 20%

# Extract the window het information

Dataset: `all50` only the finwhale individuals

```bash
# Sun Dec  6 15:05:40 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/window_het

qsub -t 1-96 step1_generate_window_het_20201206.sh

qacct -j 5655043 | grep 'exit_status' | grep -c 'exit_status  0'
```

Dataset: `f50b4` 50 finwhale indivdiduals and 4 other baleen whales

```bash
# Wed Jan 27 15:20:13 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/window_het

qsub -t 1-96 step1_generate_window_het_f50b4_20210127.sh

qacct -j 6327651 | grep 'exit_status' | grep -c 'exit_status  0'
```

Output file:

`data/window_het/raw_data/all50_window_het_20201206`

`data/window_het/raw_data/f50b4_window_het_20210128`


# Generate plotting objects

```bash
Rscript --vanilla step2_generate_plot_object_all50_20210824.R
Rscript --vanilla step2_generate_plot_object_f50b4_20210824.R
```

# Generate plots

```bash
Rscript --vanilla step3_FigS4S5_all50_window_het_20210824.R
Rscript --vanilla step3_FigS6_f50b4_winhet_barplot_20210626.R
```

# Generate summary statistics and histogram data of windowed heterozygosity

```bash
Rscript --vanilla step4_window_het_histdt_20210824.R
```





