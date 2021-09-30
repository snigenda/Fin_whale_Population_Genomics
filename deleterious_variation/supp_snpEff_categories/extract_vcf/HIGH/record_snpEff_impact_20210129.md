# step1: Extract the infomation of Annotation_Impact using the snpEff annotations

```bash
# Fri Dec 18 13:29:49 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -t 1-96 snpEff_impact/step1_snpEff_impact_extract_filteredvcf_bed_20210129.sh
qacct -j 6358785 | grep 'exit_status' | grep -c 'exit_status  0'
```

# step2: Merge the bedfiles

```bash
# Sat Jan 30 00:11:42 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub snpEff_impact/step2_snpEff_impact_merge_filteredvcf_bed_20210129.sh
qacct -j 6362230

# # Fix bugs, delete files
# rm -r ${WORKDIR}
# cd ${SCRATCHDIR}
# rm *filteredvcf_snpEff_sorted.bed
```

Output statistics:

```
[2021-01-30 00:50:33] Getting total length of JointCalls_all50_filterpass_MODIFIER_filteredvcf_snpEff.bed
20007246
[2021-01-30 00:50:57] Getting total length of JointCalls_all50_filterpass_LOW_filteredvcf_snpEff.bed
174197
[2021-01-30 00:50:58] Getting total length of JointCalls_all50_filterpass_MODERATE_filteredvcf_snpEff.bed
124824
[2021-01-30 00:50:58] Getting total length of JointCalls_all50_filterpass_HIGH_filteredvcf_snpEff.bed
2993
```

# step3: extract the filteredvcf

```bash
# Sun Jan 31 19:01:19 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub snpEff_impact/step3_snpEff_impact_extract_filteredvcf_vcf_20210129.sh "HIGH" # JOB ID: 6391119
qsub snpEff_impact/step3_snpEff_impact_extract_filteredvcf_vcf_20210129.sh "MODERATE" # JOB ID: 6391738
qsub snpEff_impact/step3_snpEff_impact_extract_filteredvcf_vcf_20210129.sh "LOW" # JOB ID: 6391744
qsub -l highp,h_data=15G,h_vmem=16G,h_rt=72:00:00 snpEff_impact/step3_snpEff_impact_extract_filteredvcf_vcf_20210129.sh "MODIFIER" # JOB ID: 6391734
```

# step3.x1: change the strategy for the modifier type

Run within individual sites

```bash
# Fri Mar 26 13:47:22 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 snpEff_impact/step3.x1_snpEff_impact_extract_filteredvcf_vcf_MODIFIER_20210326.sh
# JOB ID: 7097790
```

Still takes a long time, up the memory requests

```bash
# Sun Mar 28 20:24:38 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 snpEff_impact/step3.x1_snpEff_impact_extract_filteredvcf_vcf_MODIFIER_20210326.sh
# JOB ID: 7099042

grep 'JOB ID 7099042' /u/project/rwayne/meixilin/fin_whale/analyses/reports/snpEff_impact/step3.x1_snpEff_impact_extract_filteredvcf_vcf_MODIFIER_20210326.out.txt | grep 'Done' | wc -l

# SUCCESS!
```

Transfer the files back to project rwayne

```bash
# Thu Apr 15 11:43:11 2021
cd /u/scratch/m/meixilin/finwhale/analyses/snpEff_impact/all50/Minke
PROJECT=/u/project/rwayne/meixilin/fin_whale/analyses/snpEff_impact/all50/Minke/MODIFIER_vcf/
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/snpEff_impact/all50/Minke/

rsync -ahv ${SCRATCHDIR} ${PROJECT}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
fi
```

# extract within the CDS region

This section is reorganized into `get_ALLregions_CDS/snpEff_impactCDS`, which used the

JointCalls_all50_08_B_VariantFiltration_ALLregions_all.vcf.gz file to be more consistent with the SIFT annotation types.

# transfer files

```bash
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/snpEff_impact/all50/Minke/

rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}
```

# step4: generate vcfR objects

```bash
# Fri Apr 16 16:19:00 2021
run_step4() {
local VCFPREFIX=${1}
local VCFFILE=${2}
local TODAY=$(date "+%Y%m%d")
local LOG=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke/logs/step4_createvcfR_${VCFPREFIX}_all50_Minke_filteredvcf_${TODAY}.log
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/snpEff_impact/step4_createvcfR_20210416.R ${VCFPREFIX} ${VCFFILE} &> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done"
}

set -o pipefail
VCFDIR=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke
run_step4 'HIGH' ${VCFDIR}/JointCalls_all50_filterpass_HIGH_filteredvcf_snpEff.vcf.gz
run_step4 'MODERATE' ${VCFDIR}/JointCalls_all50_filterpass_MODERATE_filteredvcf_snpEff.vcf.gz
run_step4 'LOW' ${VCFDIR}/JointCalls_all50_filterpass_LOW_filteredvcf_snpEff.vcf.gz
```

## sync back the files

```bash
# Fri Apr 16 16:35:40 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke/
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/snpEff_impact/all50/Minke

rsync -ahv --update --exclude ".DS_Store" --exclude "Icon*" -e "ssh -i ~/.ssh/id_hoff_rsa" ${LOCAL} ${REMOTE}
```

## generate vcfR files for modifier vcfs

### first install vcfR in hoffman2

```R
# Fri Apr 16 18:54:27 2021
.libPaths(new='/u/home/m/meixilin/R/x86_64-pc-linux-gnu-library/3.6')
dir.create('/u/home/m/meixilin/R/x86_64-pc-linux-gnu-library/3.6', showWarnings = FALSE, recursive = TRUE)
install.packages('vcfR')
```
### use hoffman2 to generate vcfR objects

```bash
# Mon Apr 19 13:10:01 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -N snpEff_impact_vcfR_MODIFIER -t 1-96 snpEff_impact/step4.x1_snpEff_impact_createvcfR_MODIFIER_20210416.sh

# JOB ID: 7371194

# Mon Apr 19 17:50:41 2021
qacct -j 7371194 | grep 'exit_status' | grep -c 'exit_status  0'
```

The md5sum of local MODIFIER and remote MODIFIER vcfR was not the same but the files were identical.

```R
dt1 = readRDS(<remotefile>)
dt2 = readRDS(<localfile>)
identical(dt1,dt2)
# TRUE
```

### merge the vcfR objects

```bash
# Mon Apr 19 21:30:22 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub snpEff_impact/step4.x2_snpEff_impact_mergevcfR_MODIFIER_20210419.sh

# JOB ID: 7375640
```

## sync the files to local

```bash
# Mon Apr 19 21:48:25 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/snpEff_impact/all50/Minke/

rsync -ahv --update --exclude "derive_data/vcfR/MODIFIER_vcf/" --exclude "derive_data/gttable/MODIFIER_vcf" -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}
```


******

**This section onwards is archived**

Did not take into consideration for all alternative alleles

```bash
# Thu Apr 15 11:47:09 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/snpEff_impact/all50/Minke

mv count_sites_20210201 archive/
```

# step4: generate count sites

```bash
# Mon Feb  1 02:33:44 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub snpEff_impact/step4_snpEff_impact_countSites_20210201.sh "HIGH"
qsub snpEff_impact/step4_snpEff_impact_countSites_20210201.sh "MODERATE"
qsub snpEff_impact/step4_snpEff_impact_countSites_20210201.sh "LOW"

```

# Transfer files

```bash
# Mon Feb  1 02:40:34 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke/
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/snpEff_impact/all50/Minke/count_sites_20210201

rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}
```


```bash
# Fri Apr 16 17:07:05 2021
run_step4() {
local VCFPREFIX=${1}
local VCFFILE=${2}
local TODAY=$(date "+%Y%m%d")
local LOG=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke/logs/step4_createvcfR_${VCFPREFIX}_all50_Minke_filteredvcf_${TODAY}.log
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/snpEff_impact/step4_createvcfR_20210416.R ${VCFPREFIX} ${VCFFILE} &> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done"
}

set -o pipefail
VCFDIR=/Users/linmeixi/google_drive/finwhale/analyses/snpEff_impact/all50/Minke

for ii in {01..02}; do
IDX=$(printf %02d ${ii})
echo ${IDX}
run_step4 "MODIFIER_${IDX}" "${VCFDIR}/MODIFIER_vcf/JointCalls_all50_filterpass_MODIFIER_filteredvcf_snpEff_${IDX}.vcf.gz"
done
```

