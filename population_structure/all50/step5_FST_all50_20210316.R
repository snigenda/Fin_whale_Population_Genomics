# Title: Analyze Fst using the new LDPruned and output tables
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Fri Mar 19 00:42:32 2021
# execute with:
# source(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/all50/step5_FST_all50_20210316.R', echo = TRUE, max.deparse.length = 1000)

# preparation --------
rm(list = ls())
cat("\014")

library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)

source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R')

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

setwd('/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/all50')
outdir = './derive_data/'
plotdir = './plots/'
logdir = './logs/'

dir.create(outdir)
dir.create(plotdir)
dir.create(logdir)

sink(file = paste0(logdir, 'FST_all50_Minke_', today, '.log'))
dataset = 'all50'
ref = 'Minke'
mafcut = '10'
gdsfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds'

sessionInfo()
getwd()

# def functions --------
get_allFST <- function(genofile, pop_code) {
    fst_pop <- snpgdsFst(gdsobj = genofile, population=as.factor(pop_code),
                         method="W&C84", autosome.only = F)
    print(paste0('Fst: ', fst_pop$Fst))
    print(paste0('Mean Fst: ',fst_pop$MeanFst))
    return(fst_pop)
}

get_pairFST <- function(genofile, pop.map, popcolname, poporder) {
    pop_code = pop.map[,popcolname]
    fst <- data.frame(t(combn(poporder,2)), stringsAsFactors = FALSE)
    colnames(fst) <- c("Loc1","Loc2")
    fst$Fst <- NA
    fst$MeanFst <- NA
    for(ii in 1:nrow(fst)){
        pair <-  as.character(fst[ii,c("Loc1","Loc2")])
        popid <- which(pop.map[,popcolname] %in% pair)
        sec <- pop.map[popid, c('SampleId', popcolname)]
        colnames(sec) <- c('sample', 'location')
        fstR <- snpgdsFst(gdsobj = genofile, sample.id = sec$sample,population=as.factor(sec$location), method="W&C84",  autosome.only = F, verbose = F, remove.monosnp = T,missing.rate = 0.2)
        fst[ii,3] <- fstR$Fst
        fst[ii,4] <- fstR$MeanFst
    }
    fstc = reshape2::dcast(data = fst, formula = Loc1 ~ Loc2, value.var = 'Fst') %>%
        tibble::column_to_rownames(var = 'Loc1')
    fstc = fstc[poporder[-length(poporder)],poporder[-1]]
    print(round(fstc, digits = 4))
    output = list(fst, fstc)
    return(output)
}

# load data --------
genofile = SNPRelate::snpgdsOpen(filename = paste0(ref, '/', gdsfile))

#List of sample ids
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# sample designations
popmap = read.csv(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/popmap_all50.csv', stringsAsFactors = FALSE) %>%
    dplyr::mutate(SubPopId2 = ifelse(SubPopId %in% c('WA', 'OR', 'BC'), 'MENP', SubPopId))
pop_code <- popmap[,'PopId']
subpop_code <- popmap[,'SubPopId']
subpop2_code <- popmap[, 'SubPopId2']

# main --------
# between groups ========
fst_pop <- get_allFST(genofile = genofile, pop_code = pop_code)
fst_subpop <- get_allFST(genofile = genofile, pop_code = subpop_code)
fst_subpop2 <- get_allFST(genofile = genofile, pop_code = subpop2_code)

# pairwise Fst ========
fstpair_subpop <- get_pairFST(genofile, popmap, popcolname = 'SubPopId', poporder = c("AK", "BC", "WA", "OR", "CA", "GOC"))
fstpair_subpop2 <- get_pairFST(genofile, popmap, popcolname = 'SubPopId2', poporder = c("AK", "MENP", "CA", "GOC"))

# output files --------
write.csv(x = fstpair_subpop[[1]], file = paste0(outdir,'pairFST_subpop_maf', mafcut, '_all50_Minke_', today, '.csv'), row.names = FALSE)
write.csv(x = fstpair_subpop[[2]], file = paste0(outdir,'pairFST_subpopM_maf', mafcut, '_all50_Minke_', today, '.csv'), row.names = FALSE)

write.csv(x = fstpair_subpop2[[1]], file = paste0(outdir,'pairFST_subpop2_maf', mafcut, '_all50_Minke_', today, '.csv'), row.names = FALSE)
write.csv(x = fstpair_subpop2[[2]], file = paste0(outdir,'pairFST_subpop2M_maf', mafcut, '_all50_Minke_', today, '.csv'), row.names = FALSE)

# cleanup --------
closefn.gds(genofile)
date()
sink()
closeAllConnections()
