# Title:Generate summary statistics for the fin whale summary
# Author: Meixi Lin
# Date: Tue Apr  5 09:54:37 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE,crayon.enabled = FALSE)

setwd('/Volumes/GoogleDrive-107286704221868532902/My Drive/finwhale/analyses/revisions_SLiM_final/')
sink(file = './derived_data/logs/step2_generate_summary_stats_maintext_20220405.log')

library(dplyr)

# def functions --------
get_mean2pop <- function(plotdt, measurevar) {
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

get_median2pop <- function(plotdt, measurevar) {
    mediandt = plotdt %>%
        dplyr::group_by(popmodel) %>%
        dplyr::summarise(medianval = median(eval(as.name(measurevar))),
                         nmeasure = n(),
                         .groups = 'drop')
    # percent change based on enp
    enpmedian = mediandt$medianval[mediandt$popmodel == 'ENP']
    mediandt = mediandt %>%
        dplyr::mutate(enpchange = (medianval - enpmedian)/enpmedian)
    # percent change based on goc migration scenario
    gocmigmedian = mediandt$medianval[mediandt$popmodel == 'GOC w/ Mig']
    mediandt = mediandt %>%
        dplyr::mutate(gocmchange = (medianval - gocmigmedian)/gocmigmedian)
    return(mediandt)
}

# mean and median in ENP alone ========
get_meanENP <- function(plotdt, measurevar) {
    meandt = plotdt %>%
        dplyr::group_by(model) %>%
        dplyr::summarise(meanval = mean(eval(as.name(measurevar))),
                         nmeasure = n(),
                         .groups = 'drop')
    # percent change based on enp before bottleneck
    enpmean = meandt$meanval[meandt$model == 'pre-bott']
    meandt = meandt %>%
        dplyr::mutate(enpchange = (meanval - enpmean)/enpmean)
    return(meandt)
}

get_medianENP <- function(plotdt, measurevar) {
    mediandt = plotdt %>%
        dplyr::group_by(model) %>%
        dplyr::summarise(medianval = median(eval(as.name(measurevar))),
                         nmeasure = n(),
                         .groups = 'drop')
    # percent change based on enp before bottleneck
    enpmedian = mediandt$medianval[mediandt$model == 'pre-bott']
    mediandt = mediandt %>%
        dplyr::mutate(enpchange = (medianval - enpmedian)/enpmedian)
    return(mediandt)
}

# def variables --------
# variables to collect
vars = c('meanHet','FROH_1Mb', 'geneticLoad', 'avgStrDel', 'avgModDel', 'avgWkDel')

sessionInfo()

# load data --------
plot2popdt = readRDS(file = './derived_data/plotdt/slim_plotdt_2pop_ancestralChange_20220405.rds')
plotenpdt = readRDS(file = './derived_data/plotdt/slim_plotdt_ENP_20220405.rds')

# main --------
# group and output 2popdt ========
for (ii in vars) {
    print(ii)
    print(get_mean2pop(plot2popdt,ii))
    print(get_median2pop(plot2popdt,ii))
}

# group and output enpdt ========
for (ii in vars) {
    print(ii)
    print(get_meanENP(plotenpdt,ii))
    print(get_medianENP(plotenpdt,ii))
}

# output files --------

# cleanup --------
date()
sink()
closeAllConnections()
