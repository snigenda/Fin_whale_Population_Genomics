# Title: Summarize final slim output files
# Author: Meixi Lin
# Date: Tue Apr  5 09:32:26 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE,stringsAsFactors = FALSE)

setwd('/Volumes/GoogleDrive-107286704221868532902/My Drive/finwhale/analyses/revisions_SLiM_final/')

library(dplyr)
library(stringr)
library(ggplot2)


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
get_lastgen <- function(slimdt, lastgen) {
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
    P1cols = str_detect(colnames(enpdt),'P1')
    colnames(enpdt)[P1cols] = str_replace(colnames(enpdt)[P1cols], pattern = 'P1', replacement = '')
    gocdt = lastdt %>% 
        dplyr::select(rep,model,ends_with('P2')) %>%
        dplyr::mutate(pop = 'GOC') %>%
        dplyr::relocate(pop, .after = 'model')
    P2cols = str_detect(colnames(gocdt),'P2')
    colnames(gocdt)[P2cols] = str_replace(colnames(gocdt)[P2cols], pattern = 'P2', replacement = '')
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
outdir = './derived_data/'
plotdir = './plots/'
dir.create(outdir)
dir.create(plotdir)

nrep = 25
finalmodels = c('finwhale_2pop_AncChange_AsymMig_20220316', 'finwhale_2pop_AncChange_NoMig_20220316',
                'finwhale_ENP_3Epoch_20220316','finwhale_ENP_3Epoch_recovery_20220316')

# list of final files
migfiles = paste0('finwhale_2pop_AncChange_AsymMig_20220316/rep', 1:nrep, 
                  '/slimfinal_finwhale_2pop_AncChange_AsymMig_20220316_rep',1:nrep,'.txt')
nomigfiles = paste0('finwhale_2pop_AncChange_NoMig_20220316/rep', 1:nrep, 
                    '/slimfinal_finwhale_2pop_AncChange_NoMig_20220316_rep',1:nrep,'.txt')
enpfiles = paste0('finwhale_ENP_3Epoch_20220316/rep', 1:nrep, 
                  '/slimfinal_finwhale_ENP_3Epoch_20220316_rep',1:nrep,'.txt')
enprecfiles = paste0('finwhale_ENP_3Epoch_recovery_20220316/rep', 1:nrep, 
                     '/slimfinal_finwhale_ENP_3Epoch_recovery_20220316_rep',1:nrep,'.txt')

# variables to collect
plotting_vars = c('meanHet', 'FROH_1Mb', 'geneticLoad',
                  'avgStrDel','avgModDel','avgWkDel')

# load data --------
migdt = read_slimreps(migfiles, nrep = nrep) %>%
    dplyr::mutate(model = 'w/ Mig')
nomigdt = read_slimreps(nomigfiles, nrep = nrep) %>%
    dplyr::mutate(model = 'w/o Mig')
enpdt = read_slimreps(enpfiles, nrep = nrep) 
enprecdt = read_slimreps(enprecfiles, nrep = nrep)

# last generation in the migration scenarios
lastmigdt = get_lastgen(migdt, lastgen = 167822)
lastnomigdt = get_lastgen(nomigdt, lastgen = 167822)

# enp generations after bottleneck
enp_prebot = get_lastgen(enpdt, lastgen = 169422) %>% dplyr::mutate(model = 'pre-bott')
enp_2gen = get_lastgen(enpdt, lastgen = 169424) %>% dplyr::mutate(model = '2 gens')
enp_20gen = get_lastgen(enpdt, lastgen = 169442) %>% dplyr::mutate(model = '20 gens')
enp_20genrec = get_lastgen(enprecdt, lastgen = 169442) %>% dplyr::mutate(model = '20 gens w/ recov')

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
# dir.create(outdir, 'allgen')
saveRDS(migdt, file = paste0(outdir, 'allgen/', finalmodels[1], '_allgen_20220405.rds'))
saveRDS(nomigdt, file = paste0(outdir, 'allgen/', finalmodels[2], '_allgen_20220405.rds'))
saveRDS(enpdt, file = paste0(outdir, 'allgen/', finalmodels[3], '_allgen_20220405.rds'))
saveRDS(enprecdt, file = paste0(outdir, 'allgen/', finalmodels[4], '_allgen_20220405.rds'))

# last generations
# dir.create(outdir, 'lastgen')
saveRDS(last2popdt, file = paste0(outdir, 'lastgen/', 'finwhale_2pop_AncChange_AsymMig_20220316_lastgen_20220405.rds'))
saveRDS(lastenpdt, file = paste0(outdir, 'lastgen/', 'finwhale_ENP_3Epoch_20220316_lastgen_20220405.rds'))
write.csv(last2popdt, file = paste0(outdir, 'lastgen/', 'finwhale_2pop_AncChange_AsymMig_20220316_lastgen_20220405.csv'))
write.csv(lastenpdt, file = paste0(outdir, 'lastgen/', 'finwhale_ENP_3Epoch_20220316_lastgen_20220405.csv'))

# plotting data frame
# dir.create(outdir, 'plotdt')
saveRDS(plot2popdt, file = paste0(outdir, 'plotdt/', 'slim_plotdt_2pop_ancestralChange_20220405.rds'))
saveRDS(plotenpdt, file = paste0(outdir, 'plotdt/', 'slim_plotdt_ENP_20220405.rds'))
write.csv(plot2popdt, file = paste0(outdir, 'plotdt/', 'slim_plotdt_2pop_ancestralChange_20220405.csv'))
write.csv(plotenpdt, file = paste0(outdir, 'plotdt/', 'slim_plotdt_ENP_20220405.csv'))


save.image(file = paste0(outdir, 'SLiM_final_step1_summary_output_20220405.RData'))

# cleanup --------

date()
closeAllConnections()



