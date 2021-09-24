Author: Meixi Lin
Last updated: Tue Aug 24 14:50:10 2021

# Qualimap summary of the mapped bam files

```bash
bash preprocessing_qualimap2csv.sh
bash preprocessing_qualimap2csv_baleen_20210321.sh
```

Output file:

`data/genome_stats/raw_data/all50_qualimap_summary_20210824.txt`
`data/genome_stats/raw_data/f50b4_qualimap_summary_20210321.txt`

# Gather information from filtered vcf files

Dataset: `all50`

```bash
# Sun Nov 15 15:36:02 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/Summary_stats/count_sites

qsub -t 96 wrapper_countSitesPerIndividual_20201115.sh
qsub -t 1-95 wrapper_countSitesPerIndividual_20201115.sh

qacct -j 5296985 | grep 'exit_status' | grep -c 'exit_status  0'
```

Dataset: `f50b4`

```bash
# Wed Jan 27 13:43:04 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 Summary_stats/count_sites/wrapper_f50b4_countSitesPerIndividual_20210127.sh

qacct -j 6313694 | grep 'exit_status' | grep -c 'exit_status  0'
```

Output file:

`data/genome_stats/raw_data/all50_count_sites_20201115`
`data/genome_stats/raw_data/f50b4_count_sites_20210127`


# Generate genomewide heterozygosity info

```bash
Rscript --vanilla plot_count_sites_20201115.R
Rscript --vanilla plot_count_sites_f50b4_20210128.R
```
