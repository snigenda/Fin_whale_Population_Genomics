# step1: Extract the infomation of LOF using the snpEff annotations

```bash
# Fri Dec 18 13:29:49 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub -t 1-96 get_ALLregions_CDS/snpEff_LOF/step1_snpEff_LOF_extract_ALLregions_bed_20201218.sh
qacct -j 5885343 | grep 'exit_status' | grep -c 'exit_status  0'
```

# step2: merge the bedfiles
```bash
# Mon Dec 28 15:25:38 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub get_ALLregions_CDS/snpEff_LOF/step2_snpEff_LOF_merge_ALLregions_bed_20201218.sh
qacct -j 6000745

# Total length: 2237 by the ALLregions
```


# step3: Extract the vcf
```bash
# Mon Dec 28 17:51:11 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub get_ALLregions_CDS/snpEff_LOF/step3_snpEff_LOF_extract_ALLregions_vcf_20201218.sh
qacct -j 6001753

# Output: /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.vcf.gz
# Total after trimAlternatives: 2164
# The other 73 sites not included in the vcf but in the bedfile did not contain variant sites after filtering, so they got removed
```

# check what was discarded in the vcf compared with the bedfiles

```bash
# Mon Dec 28 20:44:10 2020
cd ~/flashscratch/playground

bedtools intersect -a /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/bedfiles/JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.bed -b /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.vcf.gz -v

gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
--preserveAlleles \
-V /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/JointCalls_all50_08_B_VariantFiltration_ALLregions_all.vcf.gz \
-o test.vcf.gz \
-L test.bed

bcftools query -f '%FILTER\n' test.vcf.gz
/u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke_canonical_CDS/subset_ALLregions_vcf_20201205/JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.vcf.gz | sort | uniq -c
bcftools stats -s- test.vcf.gz
```

# Transfer the file locally

```bash
# Tue Feb  2 11:40:16 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/get_ALLregions_CDS/all50/Minke/
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/JointCalls_all50_filterpassmiss*.vcf.gz*

rsync -ahv --checksum -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}
```

# Get LOF types

|  | LOFtype | Count |
| :--- | :---------------------------------------------------------- | :--- |
| 1    | splice_acceptor_variant&intron_variant                      | 376  |
| 2    | splice_acceptor_variant&splice_donor_variant&intron_variant | 73   |
| 3    | splice_donor_variant&intron_variant                         | 547  |
| 4    | start_lost                                                  | 212  |
| 5    | start_lost&splice_region_variant                            | 2    |
| 6    | stop_gained                                                 | 923  |
| 7    | stop_gained&splice_region_variant                           | 31   |

# Difference with the previous version:

Added the splice region variants.

# Conclusions:
1. The tally of two LOF types does not change the relative percentage and performance between sites
2. Some LOF does not fall within the canonical coding regions but are STILL very important.
3. Some sites are included in the bedfile but did not end up in the vcf file because the site was not a SNP anymore after WGSproc8 filtering



