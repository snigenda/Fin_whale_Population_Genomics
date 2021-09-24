# step1: output CDS

```bash
# Fri Dec 18 09:35:44 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub get_ALLregions_CDS/gtf_CDS/step1_output_CDS_from_gtf_20201217.sh

qacct -j 5883119
```

# step2: extract vcf from these regions

```bash
# Fri Dec 18 09:38:02 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub get_ALLregions_CDS/gtf_CDS/step2_extract_ALLregions_CDS_vcf_20201218.sh

qacct -j 5883171
```

* bcftools stats

```
[2020-12-18 15:02:23] Concatenating vcf to /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/JointCalls_all50_08_B_VariantFiltration_ALLregions_all.vcf.gz ...
Number of samples: 50
Number of SNPs:    528923
Number of INDELs:  27730
Number of MNPs:    0
Number of others:  0
Number of sites:   30896030
```

# check that the regions match in the vcf and bedfiles

```bash
check_length() {
    local FILENAME=${1}
    echo -e "[$(date "+%Y-%m-%d %T")] Getting total length of ${FILENAME}"
    awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${FILENAME}
}

cd /u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke
# subset the bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.bed to the regions with Haplotypes called
Rscript --vanilla /u/project/rwayne/meixilin/fin_whale/analyses/scripts/get_ALLregions_CDS/gtf_CDS/x1_subset_bed2contiglist_20201218.R
echo $?

check_length bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.inVCF.bed
# 30896030

# check that the files were matching
bedtools intersect -sorted -a bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.inVCF.bed -b JointCalls_all50_08_B_VariantFiltration_ALLregions_all.vcf.gz -v
bedtools intersect -sorted -a JointCalls_all50_08_B_VariantFiltration_ALLregions_all.vcf.gz -b bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.inVCF.bed -v
# YES!
```