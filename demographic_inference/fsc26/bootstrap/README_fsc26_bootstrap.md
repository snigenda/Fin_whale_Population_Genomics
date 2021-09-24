[TOC]

# step0: create bootstrap files 

Modifies the `*_maxL.par` file in the result folder to generate DNA sequence data for bootstrap estimation. 

* Use the best likelihood estimates from previous runs in the  [fsc_output_list.csv](fsc_output_list.csv). The runs marked with `#` were not included in the final SI table.

```bash
bash step0_generate_btsp_files_20210603.sh 'local' 'neutral' '15'
```

Output: This step generates the files listed in the [neutral](neutral) folder. For example, for the `15` (2D.SplitAsyMigTw2) model, you obtain these four files

```
2D.SplitAsyMigTw2.GOCENP_vreal_r20_boot.est 
2D.SplitAsyMigTw2.GOCENP_vreal_r20_boot.tpl
2D.SplitAsyMigTw2.GOCENP_vreal_r20_boot.pv
2D.SplitAsyMigTw2.GOCENP_vreal_r20_boot.par
```

# step1: generate bootstrap SFS

Using the `*_boot.par` file generated, you create 100 bootstrap SFS in the given demographic scenario. 

```bash
# generate 100 bootstrap SFS for 2D.SplitAsyMigTw2 model
qsub -N btspSFS_15 step1_bootstrapSFS_20210603.sh '15'
```

# step2: run fastsimcoal estimations

Run fastsimcoal estimates on the bootstrap SFS. 

```bash
# For one population (1D) models
qsub step2_runfsc_btSFS_1DModels_n20rd_t200_20210705.sh ${REFID}
# For two population (2D) models
qsub step2_runfsc_btSFS_2DModels_n20rd_20210823.sh ${REFID}
```



# step3: gather results from each runs

This step parses the best likelihood estimates from the 20 replicates of each of the 100 bootstrap SFSes and add the best likelihood estimates from the data as line0. Generates a csv table with 2000 bootstrap + 1 real-life data observations. 

```
${WORKDIR}/boot_${ii}/run_${jj}/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot.bestlhoods
```

# step4: format the results summary used in the Supplemental Information


```bash
# Tue Sep  7 23:58:32 2021
LOCAL=~/Lab/finwhale_manuscript/data/demography/fsc26/bootstrap/raw_data/
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/fsc26/param_btsp/neutral/resultsSummaries/./*/n20rd_202109*

rsync -ahv --update --relative -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}

cd ${LOCAL}
wc -l */*/*_Summary_*.csv
```

Use the R script to convert values to diploids and get intervals. 

```R
source('/Users/linmeixi/Lab/finwhale_manuscript/scripts/demography/fsc26/bootstrap/fsc_summary_results.R')
```

