# step0: generate SNPRelate files

```bash
# Tue Mar 16 11:09:08 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub PopStructure/all50/wrapper_step0_vcf2gds_all50_20210316.sh

qacct -j 6885061

```

# step1: gds filter for PASS sites and LD pruning

```bash
# Wed Mar 17 11:29:42 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub PopStructure/all50/wrapper_step1_gdsLDPruning_all50_20210316.sh

qacct -j 6908298

cd /u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/all50/Minke
gzip JointCalls_all50_filterpass_bialleic_all_20210317.txt
```

# step2: gds file to vcf formats

```bash
# Thu Mar 18 14:29:58 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub PopStructure/all50/wrapper_step2_Prunedgds2vcf_all50_20210316.sh

qacct -j 6933390
```

# step3: convert vcf to plink format

```bash
# Thu Mar 18 15:18:26 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/all50/Minke
WORKSCRIPT=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/PopStructure/all50/step3_Prunedvcf2plink_all50_20210316.sh
bash ${WORKSCRIPT} 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf05_SA_mrF.vcf' 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf05_SA_mrF'
bash ${WORKSCRIPT} 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10_SA_mrF.vcf' 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10_SA_mrF'
bash ${WORKSCRIPT} 'JointCalls_all50_filterpass_bialleic_all_LDPruned_mafNA_SA_mrF.vcf' 'JointCalls_all50_filterpass_bialleic_all_LDPruned_mafNA_SA_mrF'
```

# step4: run admixture

Run all the maf cutoffs

```bash
# Thu Mar 18 16:11:40 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub PopStructure/all50/step4_admixture_all50_20210316.sh

qacct -j 6937177
```

# stepx: transfer the files locally

```bash
# Thu Mar 18 16:12:29 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/all50/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/all50/Minke/

rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" --exclude 'Admixture_20210318' ${REMOTE} ${LOCAL}

# Fri Mar 19 00:02:41 2021
rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" --delete --update ${REMOTE} ${LOCAL}

```
**IMPORTANT NOTE: From step4.1 here on the maf cutoff is set at 0.10**

```R
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
mafcut = '10'
Klist = 2:6
subpoporder = c("AK", "BC", "WA", "OR", "CA", "GOC") # order for subpopulations 

gdsfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds'
# confirm the sample names
plinkfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10_SA_mrF.nosex'
```

# step4.1: plotting admixture

# step6: run relatedness analyses

```R
# Tue Apr  6 11:18:26 2021
ENPibd.coeff = read.csv(file = 'derive_data/ENPibdcoeff_maf10_SUBmaf05miss05_all50_Minke_20210323.csv', stringsAsFactors = FALSE)
summary(ENPibd.coeff$kinship)

GOCibd.coeff = read.csv(file = 'derive_data/GOCibdcoeff_maf10_SUBmaf05miss05_all50_Minke_20210323.csv', stringsAsFactors = FALSE)
summary(GOCibd.coeff$kinship)

# Also ran wilcoxon test
```
