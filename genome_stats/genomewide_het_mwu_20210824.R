# Title: Mann Whitney U test for the genomewide heterozygosity
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Aug 24 20:33:02 2021

# Relevant values: In GOC individuals we found patterns of reduced variation, with an average 1.13 heterozygotes per kb ... In contrast, the ENP population had much higher diversity (1.76 het/kb; Mann-Whitney U [MWU] test p<0.001; Figure 2A) 

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
library(dplyr)
library(reshape2)

# def functions --------
# pull out histogram data


# def variables --------
setwd('~/Lab/finwhale_manuscript/')

# load data --------
totaldt = read.csv(file = './data/genome_stats/derived_data/all50_genomewide_heterozygosity_20210824.csv', row.names = 1, stringsAsFactors = FALSE)

# main --------
# mwu test
wilcox.test(GenomeHet ~ PopId, data = totaldt)
# Wilcoxon rank sum test
# 
# data:  GenomeHet by PopId
# W = 580, p-value = 1.152e-10
# alternative hypothesis: true location shift is not equal to 0

# mean values
meandt = totaldt %>%
    dplyr::group_by(PopId) %>%
    dplyr::summarise(PopGenomeHet = mean(GenomeHet),
                     .groups = 'drop')
meandt
# A tibble: 2 x 2
#   PopId PopGenomeHet
#   <chr>        <dbl>
# 1 ENP        0.00176
# 2 GOC        0.00113

# output files --------

# cleanup --------
date()
closeAllConnections()
