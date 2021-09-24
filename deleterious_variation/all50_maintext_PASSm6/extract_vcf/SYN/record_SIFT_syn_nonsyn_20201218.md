# step1: extract the synonymous/nonsynonymous regions

```bash
# Fri Dec 18 09:51:13 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -t 1-96 get_ALLregions_CDS/SIFT_syn_nonsyn/step1_SIFT_syn_nonsyn_extract_ALLregions_bed_20201218.sh

qacct -j 5883178 | grep 'exit_status' | grep -c 'exit_status  0'
```

# step2: merge the regions

```bash
# Fri Dec 18 11:10:45 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub get_ALLregions_CDS/SIFT_syn_nonsyn/step2_SIFT_syn_nonsyn_merge_ALLregions_bed_20201218.sh

qacct -j 5884334
```

## Output:

```
[2020-12-18 11:24:46] Getting total length of /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS//all50/Minke/bedfiles/JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.bed
153989

[2020-12-18 11:24:47] Getting total length of /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS//all50/Minke/bedfiles/JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.bed
133934
```

## check that the syn/nonsyn regions all fall in the gtf files

### NOTES from Fri Aug 27 14:09:04 2021

I did extract all the bedfiles from the original filteredvcf just like what you would expect. The only difference from lof05 was that I did not pull out the filteredvcf the same way as LOF05 since I had extracted a CDS only file before.

No need to rerun the extraction. Though the extraction process could have been better optimized.

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/bedfiles
# A should completely fall in B
bedtools intersect -a JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.bed -b Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.inVCF.bed -v
bedtools intersect -a JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.bed -b Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.inVCF.bed -v
```

# step3: extract the syn/nonsyn regions

```bash
# Sat Dec 19 00:46:39 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub get_ALLregions_CDS/SIFT_syn_nonsyn/step3_SIFT_syn_nonsyn_extract_ALLregions_vcf_20201218.sh
qacct -j 5892864
```

## Output

```
JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.vcf.gz
Number of samples: 50
Number of SNPs:    152982
Number of INDELs:  0
Number of MNPs:    0
Number of others:  0
Number of sites:   152982

JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.vcf.gz
Number of samples: 50
Number of SNPs:    132525
Number of INDELs:  0
Number of MNPs:    0
Number of others:  0
Number of sites:   132525
```
