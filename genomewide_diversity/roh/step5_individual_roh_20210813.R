# Title: Plot each individual's ROH
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Mar 28 10:57:11 2021
# Modification: Add bcftools options too
# Date: Fri Aug 13 16:48:42 2021
# Modification: Change to the final directory
# Date: Sun Sep 12 10:43:27 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

library(dplyr)
library(ggplot2)

# def functions --------
calculate_froh <- function(roh, minlen = 1e+6, totallen) {
    output = roh %>%
        dplyr::filter(length >= minlen) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(froh = sum(length)/totallen, 
                         .groups = "drop")
    return(output)
}

join_frohdtl <- function(frohdtl, frohnames) {
    frohdt = dplyr::full_join(x = frohdtl[[1]], y = frohdtl[[2]], by = 'sample')
    frohdt = dplyr::full_join(x = frohdt, y = frohdtl[[3]], by = 'sample') 
    colnames(frohdt) = frohnames
    return(frohdt)
}

# def variables --------
indir = "./data/roh/derived_data/"
outdir = "./data/roh/derived_data/"

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

# main --------
frohdtl_zoo = lapply(rohlens, calculate_froh, roh = zooroh, totallen = genomelen)
frohdtl_bcf = lapply(rohlens, calculate_froh, roh = bcfroh, totallen = genomelen)

frohdt_zoo = join_frohdtl(frohdtl = frohdtl_zoo, 
                      frohnames = c('Sample', 'F_ROH_100k_zoo', 'F_ROH_1M_zoo', 'F_ROH_5M_zoo'))

frohdt_bcf = join_frohdtl(frohdtl = frohdtl_bcf, 
                          frohnames = c('Sample', 'F_ROH_100k_bcf', 'F_ROH_1M_bcf', 'F_ROH_5M_bcf'))

froh_all = dplyr::full_join(frohdt_zoo, frohdt_bcf, by = 'Sample')

# output files --------
write.csv(froh_all, file = paste0(outdir, 'FROH_100k1M5M_all_', today, '.csv'), row.names = FALSE, quote = FALSE)

# cleanup --------
date()
closeAllConnections()
