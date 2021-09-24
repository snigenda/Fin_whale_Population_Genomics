# check the two versions of RDS files

```bash
md5 -r /Users/linmeixi/Lab/finwhale_manuscript/data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_20210824.rds
# 5217686e012ad5b68f8aa2773193b1d9
md5 -r /Users/linmeixi/google_drive/finwhale/analyses/window_het/all50/derive_data/winHet_1Mbwin_1Mbstep_20Per_allsamples_plotdt_20210325.rds
# aea267d5496862447ed487fbac4a9110
```

``` R
dt1 = readRDS('/Users/linmeixi/google_drive/finwhale/analyses/window_het/all50/derive_data/winHet_1Mbwin_1Mbstep_20Per_allsamples_plotdt_20210325.rds')
dt2 = readRDS('/Users/linmeixi/Lab/finwhale_manuscript/data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_20210824.rds')

identical(dt1,dt2)

# change the dt1 file structure
dt1$chrom = as.character(dt1$chrom)

identical(dt1,dt2)
# TRUE (so no problem overall, use the new version)
```

# matching the Figure S6 observed heterozygosity 

```R
genomehet = read.csv(file = "./data/genome_stats/derived_data/f50b4_genomewide_heterozygosity_20210824.csv", 
+                          comment.char = '#',stringsAsFactors = FALSE)
het2 = read.csv('./scripts/config/Baleen_Genomewide_Het_20210906.csv')
het2$Observed.pi[2] == genomehet$GenomeHet[1]
het2$Observed.pi[4] == genomehet$GenomeHet[2]
het2$Observed.pi[5] == genomehet$GenomeHet[33]
het2$Observed.pi[6] == genomehet$GenomeHet[54]
```

# match the Figure S6 histogram

```R
dt = read.csv('/Users/linmeixi/google_drive/finwhale/analyses/window_het/f50b4/derive_data/winHet_1Mbwin_1Mbstep_20Per_allsamples_barhistdt_20210626.csv', row.names = 1)
dt = dt[c(1:11,23:44,12:22,45:264),c(2,1,4)]
colnames(dt)[1] = 'maxhetpkb_ex'
rownames(dt) = 1:264
testdt = barplotdt %>%
    dplyr::mutate(sample = as.character(sample))
all.equal(testdt, dt) # Yes, they are the same!
testdt[,1] == dt[,1]
testdt[,2] == dt[,2]
hist(testdt[,3] - dt[,3])
```

