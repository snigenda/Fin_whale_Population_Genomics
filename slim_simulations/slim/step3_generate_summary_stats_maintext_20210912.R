# Title:Generate summary statistics for the fin whale summary
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Sep 12 13:26:35 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE,crayon.enabled = FALSE)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')
# sink(file = './data/slim/derived_data/logs/step3_generate_summary_stats_maintext_20210912.log')

library(dplyr)

# def functions --------
get_meanval <- function(plotdt, measurevar) {
    meandt = plotdt %>%
        dplyr::group_by(popmodel) %>%
        dplyr::summarise(meanval = mean(eval(as.name(measurevar))),
                         nmeasure = n(),
                         .groups = 'drop')
    # percent change based on enp 
    enpmean = meandt$meanval[meandt$popmodel == 'ENP']
    meandt = meandt %>%
        dplyr::mutate(enpchange = (meanval - enpmean)/enpmean)
    # percent change based on goc migration scenario
    gocmigmean = meandt$meanval[meandt$popmodel == 'GOC w/ Mig']
    meandt = meandt %>%
        dplyr::mutate(gocmchange = (meanval - gocmigmean)/gocmigmean)
    return(meandt)
}

# def variables --------
indir = './data/slim/derived_data/'

# variables to collect
plotting_vars = c("avgStrDel","avgModDel","avgWkDel")
ylabs = c('# str.del. alleles', '# mod.del. alleles', '# weak del. alleles')

sessionInfo()

# load data --------
plot2popdt = readRDS(file = paste0(indir, 'slim_plotdt_2pop_ancestralChange_20210919.rds'))

# main --------
get_meanval(plot2popdt, measurevar = 'meanHet')
get_meanval(plot2popdt, measurevar = 'FROH_1Mb')
get_meanval(plot2popdt, measurevar = 'geneticLoad')
get_meanval(plot2popdt, measurevar = 'avgStrDel')
get_meanval(plot2popdt, measurevar = 'avgModDel')
get_meanval(plot2popdt, measurevar = 'avgWkDel')

# output files --------

# cleanup --------
date()
closeAllConnections()
