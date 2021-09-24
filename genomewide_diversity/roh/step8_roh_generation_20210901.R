# Title: Get average ROH lengths and the generation of coalescent
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Sep  1 20:36:02 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/Users/linmeixi/google_drive/finwhale/analyses/important_results/Runs_of_homozygosity/")

library(dplyr)
library(ggplot2)

# def functions --------

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

# load data --------
zooroh = readRDS(file = './derived_data/ROH_zooroh_summary_final_20210812.rds')

# get average length
meanroh = zooroh %>%
    dplyr::group_by(rohcat, sample, pop) %>%
    dplyr::summarise(meanlen = mean(length), .groups = 'drop') 

poproh = meanroh %>%
    dplyr::group_by(rohcat, pop) %>%
    dplyr::summarise(poplen = mean(meanlen)) %>%
    dplyr::mutate(inbred_gen = 100/(2*poplen/1e+6))

# main --------

# output files --------
write.csv(poproh, file = paste0('./derived_data/meanroh_zoo_', today, '.csv'))

# cleanup --------
date()
closeAllConnections()
