# Title: Summarize final slim output files
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Sep 10 22:54:25 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE,stringsAsFactors = FALSE)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

library(dplyr)
library(stringr)
library(ggplot2)

sink(file = './data/slim/derived_data/logs/summary_output_20210910.log')

# def functions --------
read_slim <- function(filename, repnum) {
    if (file.exists(filename)) {
        dtlines = readLines(filename)
        startline = which(str_detect(dtlines, pattern = '^gen,') == TRUE)
        if (length(startline) == 0) {
            print(filename)
            return(NA)
        } else {
            dt = read.csv(filename, skip = startline-1) %>%
                dplyr::mutate(rep = repnum)
            if(colnames(dt)[1] != 'gen') {
                stop('Skipped wrong amount of lines')
            }
            # print the file name and replication number
            print(paste0('Rep', repnum, ':', filename))
        }
        return(dt)
    } else {
        return(NA)
    }
}

read_slimreps <- function(slimfiles, nrep) {
    slimdt = vector('list', length = nrep)
    for (ii in 1:nrep) {
        slimdt[[ii]] = read_slim(slimfiles[ii],repnum = ii)
    }
    slimdt = dplyr::bind_rows(slimdt[!is.na(slimdt)])
    return(slimdt)
}

# need to be raw slim data
get_lastgen <- function(slimdt, lastgen = 167732) {
    outdt = slimdt %>%
        dplyr::group_by(rep) %>%
        dplyr::filter(gen == lastgen) %>%
        dplyr::ungroup()
    return(outdt)
}

format_last2popdt <- function(lastdt) {
    enpdt = lastdt %>% 
        dplyr::select(rep,model,ends_with('P1')) %>%
        dplyr::mutate(pop = 'ENP') %>%
        dplyr::relocate(pop, .after = 'model')
    colnames(enpdt)[4:ncol(enpdt)] = str_replace(colnames(enpdt)[4:ncol(enpdt)], pattern = 'P1', replacement = '')
    gocdt = lastdt %>% 
        dplyr::select(rep,model,ends_with('P2')) %>%
        dplyr::mutate(pop = 'GOC') %>%
        dplyr::relocate(pop, .after = 'model')
    colnames(gocdt)[4:ncol(gocdt)] = str_replace(colnames(gocdt)[4:ncol(gocdt)], pattern = 'P2', replacement = '')
    outdt = dplyr::bind_rows(enpdt, gocdt) %>% 
        dplyr::mutate(popmodel = paste(pop, model))
    return(outdt)
}

format_lastenpdt <- function(lastdt, pop = 'ENP', 
                            varlevel = c('pre-bott','2 gens','20 gens', '20 gens w/ recov')) {
    lastdt$model = factor(lastdt$model, levels = varlevel)
    colnames(lastdt) = str_replace(colnames(lastdt), pattern = 'P1', replacement = '')
    outdt = lastdt %>%
        dplyr::mutate(pop = pop)
    return(outdt)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
outdir = './data/slim/derived_data/'
indir = './data/slim/raw_data/'

nrep = 25
finalmodels = c('finwhale_2pop_ancestralChange_040121', 'finwhale_2pop_ancestralChange_noMig_040121',
           'finwhale_ENP_022521','finwhale_ENP_recovery_022521')

# list of final files
migfiles = paste0(indir, 'finwhale_2pop_ancestralChange_040121/rep', 1:25, '/slimout_finwhale_2pop_ancestralChange_040121_rep',1:25,'.txt')
nomigfiles = paste0(indir, 'finwhale_2pop_ancestralChange_noMig_040121/rep', 1:25, '/slimout_finwhale_2pop_ancestralChange_noMig_040121_rep',1:25,'.txt')
enpfiles = list.files(path = paste0(indir, 'finwhale_ENP_022521'), pattern = 'tail', recursive = T, full.names = T)
enprecfiles = list.files(path = paste0(indir, 'finwhale_ENP_recovery_022521'), pattern = 'tail', recursive = T, full.names = T)

# variables to collect
plotting_vars = c('meanHet', 'FROH_1Mb', 'geneticLoad',
                  'avgStrDel','avgModDel','avgWkDel')

# load data --------
migdt = read_slimreps(migfiles, nrep = nrep) %>%
    dplyr::mutate(model = 'w/ Mig')
nomigdt = read_slimreps(nomigfiles, nrep = nrep) %>%
    dplyr::mutate(model = 'w/o Mig')
enpdt = read_slimreps(enpfiles, nrep = nrep) # note the nrep might not match .1 .2 measures, see logs for detail
enprecdt = read_slimreps(enprecfiles, nrep = nrep)

# last generation in the migration scenarios
lastmigdt = get_lastgen(migdt, lastgen = 167715)
lastnomigdt = get_lastgen(nomigdt, lastgen = 167715)

# enp generations after bottleneck
enp_prebot = get_lastgen(enpdt, lastgen = 169490) %>% dplyr::mutate(model = 'pre-bott')
enp_2gen = get_lastgen(enpdt, lastgen = 169492) %>% dplyr::mutate(model = '2 gens')
enp_20gen = get_lastgen(enpdt, lastgen = 169510) %>% dplyr::mutate(model = '20 gens')
enp_20genrec = get_lastgen(enprecdt, lastgen = 169510) %>% dplyr::mutate(model = '20 gens w/ recov')

# main --------
# combine the 2pop models
last2popdt = dplyr::bind_rows(lastmigdt, lastnomigdt)

# combine recovery
lastenpdt = dplyr::bind_rows(enp_prebot, enp_2gen, enp_20gen, enp_20genrec)

# format plotdt for 2pop model ========
# only keep one ENP for reference (since they are very similar)
plot2popdt = format_last2popdt(last2popdt) %>%
    dplyr::mutate(geneticLoad = 1 - meanFitness) %>%
    dplyr::filter(popmodel != 'ENP w/o Mig') %>%
    dplyr::mutate(popmodel = ifelse(popmodel == 'ENP w/ Mig', 'ENP', popmodel))

# format plotdt for enp model ========
plotenpdt = format_lastenpdt(lastenpdt) %>%
    dplyr::mutate(geneticLoad = 1 - meanFitness) 

# output files --------
# raw files
saveRDS(migdt, file = paste0(outdir, finalmodels[1], '_allgen_20210910.rds'))
saveRDS(nomigdt, file = paste0(outdir, finalmodels[2], '_allgen_20210910.rds'))
saveRDS(enpdt, file = paste0(outdir, finalmodels[3], '_allgen_20210910.rds'))
saveRDS(enprecdt, file = paste0(outdir, finalmodels[4], '_allgen_20210910.rds'))

# last generations
saveRDS(last2popdt, file = paste0(outdir, 'finwhale_2pop_ancestralChange_040121_lastgen_20210910.rds'))
saveRDS(lastenpdt, file = paste0(outdir, 'finwhale_ENP_022521_lastgen_20210910.rds'))

# plotting data frame
saveRDS(plot2popdt, file = paste0(outdir, 'slim_plotdt_2pop_ancestralChange_20210910.rds'))
saveRDS(plotenpdt, file = paste0(outdir, 'slim_plotdt_ENP_20210910.rds'))

# cleanup --------
sink()
date()
closeAllConnections()



