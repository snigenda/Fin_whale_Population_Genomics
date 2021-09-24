# Title: Calculate Froh and summarized length
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Sep  8 00:25:56 2021
# Related to main text: wilcoxon test and generate data frame for downstream plotting

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(ggplot2)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

# def functions --------
calculate_froh <- function(roh, minlen = 1e+6, totallen) {
    output = roh %>%
        dplyr::filter(length >= minlen) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(froh = sum(length)/totallen, 
                         .groups = "drop")
    return(output)
}

# calculate proportions
calculate_nlenrohcat <- function(roh, software) {
    output = roh %>%
        dplyr::group_by(sample, rohcat) %>%
        dplyr::summarise(countcat = n(), 
                         sumcat = sum(length),
                         .groups = "drop")
    outputn = output %>%
        reshape2::dcast(., sample ~ rohcat, value.var = "countcat")
    colnames(outputn)[-1] = paste0("N_", software, "_", colnames(outputn)[-1])
    outputn[is.na(outputn)] = 0
    outputlen = output %>%
        reshape2::dcast(., sample ~ rohcat, value.var = "sumcat")
    colnames(outputlen)[-1] = paste0("LEN_", software, "_", colnames(outputlen)[-1])
    outputlen[is.na(outputlen)] = 0
    output2 = dplyr::left_join(outputn, outputlen, by = "sample")
    # output both the long and short format (the long format will be used in maintext plotting)
    output1 = output %>%
        dplyr::mutate(software = software)
    return(list(output1, output2))
}

# def variables --------
outdir = "./data/roh/derived_data/"
dir.create(outdir)

today = format(Sys.Date(), "%Y%m%d")
genomelen = 2324429847

# new roh length and categories (August 2021)
rohlens = c(0.1,1,5)*1e+6
rohcats = c('0.1_1','1_5','5_Inf')
lenslab <- c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb')

sessionInfo()

# load data --------
bcfroh = readRDS(file = './data/roh/derived_data/ROH_bcfroh_summary_final_20210907.rds')
zooroh = readRDS(file = './data/roh/derived_data/ROH_zooroh_summary_final_20210907.rds')
genomehet = read.csv(file = './data/genome_stats/derived_data/all50_genomewide_heterozygosity_20210824.csv',
                     row.names = 1, stringsAsFactors = F) %>%
    dplyr::select(-starts_with("Total"))

# main --------
# get categorized ROH types ========
bcfroh_nlenrohcat2 = calculate_nlenrohcat(bcfroh, software = "bcf")
zooroh_nlenrohcat2 = calculate_nlenrohcat(zooroh, software = "zoo")

# the long format
forplot = dplyr::bind_rows(bcfroh_nlenrohcat2[[1]],zooroh_nlenrohcat2[[1]])
# the short format
bcfroh_nlenrohcat = bcfroh_nlenrohcat2[[2]]
zooroh_nlenrohcat = zooroh_nlenrohcat2[[2]]

# get Froh (not final, using only one measurement of >= 1Mb) ========
bcfroh_froh = calculate_froh(bcfroh, minlen = 1e+6, totallen = genomelen)
zooroh_froh = calculate_froh(zooroh, minlen = 1e+6, totallen = genomelen)
summary(zooroh_froh)

# write a summary table ========
bcfrohsum = dplyr::left_join(bcfroh_froh, bcfroh_nlenrohcat, by = "sample")
zoorohsum = dplyr::left_join(zooroh_froh, zooroh_nlenrohcat, by = "sample") 

rohsum = dplyr::left_join(genomehet,zoorohsum, by = c('SampleId' = 'sample')) %>%
    dplyr::left_join(., bcfrohsum, by = c('SampleId' = 'sample'), suffix = c("_zoo", "_bcf")) 
# get samples not called in bcftools
rohsum[which(is.na(rohsum$froh_bcf)), 'SampleId']
# [1] "ENPCA09" "ENPOR12" "GOC010" 

# get total N and length
rohsum2 = rohsum %>%
    dplyr::mutate(N_zoo_all = N_zoo_0.1_1 + N_zoo_1_5 + N_zoo_5_Inf,
                  N_bcf_all = N_bcf_0.1_1 + N_bcf_1_5 + N_bcf_5_Inf,
                  LEN_zoo_all = LEN_zoo_0.1_1 + LEN_zoo_1_5 + LEN_zoo_5_Inf,
                  LEN_bcf_all = LEN_bcf_0.1_1 + LEN_bcf_1_5 + LEN_bcf_5_Inf)

# run a wilcoxon test (in the main text) =========
# "Overall, GOC individuals contained considerably more ROH segments than ENP individuals (MWU test p<0.001)"
wilcox.test(N_zoo_all ~ PopId, data = rohsum2)
# Wilcoxon rank sum test with continuity correction
# 
# data:  N_zoo_all by PopId
# W = 30, p-value = 9.415e-08
# alternative hypothesis: true location shift is not equal to 0

# output files --------
# the rohsummary file -> will be used in step4, bcftools and zooroh compare
write.csv(rohsum2, file = paste0(outdir, "rohsummary_zoobcf_", today, ".csv"))
saveRDS(rohsum2,  file = paste0(outdir, "rohsummary_zoobcf_", today, ".rds"))
# the forplot file -> will be used in the main text plotting
saveRDS(forplot, file = paste0(outdir, "rohsummarylong_zoobcf_", today, ".rds"))

# cleanup --------
date()
closeAllConnections()
