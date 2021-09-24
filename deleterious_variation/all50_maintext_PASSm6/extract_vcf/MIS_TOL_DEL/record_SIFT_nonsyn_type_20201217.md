# step1: extract bedfiles

```bash
# Mon Dec 28 13:55:44 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub get_ALLregions_CDS/SIFT_nonsyn_type/step1_extract_NONSYNtype_ALLregions_bed_20201218.sh

# Mon Dec 14 17:03:27 2020
qacct -j 6000297
```

Notes from Fri Aug 27 15:57:31 2021:

In this step, the bedfile was extracted from the file that is already comprised of only nonsynonymous variations.

Have already double checked all the numbers matched.

`132525` nonsyn variations = 24922+8750+98750+103

# get a tally of file types

```bash
# Mon Dec 28 14:13:19 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

Rscript --vanilla get_ALLregions_CDS/SIFT_nonsyn_type/step1.2_tally_NONSYNtype_ALLregions_bed_20201218.R
echo $?
```

Output:

```
DELETERIOUS DELETERIOUS_(*WARNING!_Low_confidence) 
     24922                                   8750 
 TOLERATED                                   <NA> 
     98750                                    103 
```

# step2: extract vcf files

```bash
# Mon Dec 28 14:54:15 2020
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub get_ALLregions_CDS/SIFT_nonsyn_type/step2_extract_NONSYNtype_ALLregions_vcf_20201218.sh
qacct -j 6000577
```

# for the final files, check if there are overlapps

```bash
cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke_canonical_CDS/bedfiles

bedtools intersect -a JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.bed -b JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.bed

bedtools intersect -a JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.bed -b JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.bed

bedtools intersect -a JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.bed -b JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.bed

```

## conclusion:

> Tue Dec 15 13:02:47 2020

Some sites did overlap. Different softwares were painful.

```
(gentools) [meixilin@n7361 bedfiles]$ bedtools intersect -a JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.bed -b JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.bed
NW_006726621.1	997486	997487
NW_006727020.1	5007154	5007155
NW_006728239.1	2434004	2434005
NW_006729790.1	8548662	8548663
(gentools) [meixilin@n7361 bedfiles]$ bedtools intersect -a JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT.bed -b JointCalls_all50_filterpassmiss_LOF05_ALLregions_all_snpEff.bed
NW_006725354.1	34976368	34976369
NW_006726399.1	1855707	1855708
NW_006727020.1	9524796	9524797
NW_006728239.1	6864650	6864651
NW_006731012.1	1703214	1703215
NW_006733789.1	7586176	7586177
```

No corrections were made since the count was minute compared to the total count and the effects take place on different transcripts in these cases.
