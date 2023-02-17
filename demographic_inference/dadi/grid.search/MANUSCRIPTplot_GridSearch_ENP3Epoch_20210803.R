# Title: Plot grid search result for the manuscript
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Aug 25 23:53:14 2021
# Modification: Add source_data
# Date: Mon Jan 23 10:48:25 2023


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd('/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/dadi/grid.search/')

library(ggplot2)
library(dplyr)
# install.packages('ggforce')
require(ggforce)
require(ggpubr)
require(ggpmisc)

# def functions --------
change_labs <- function(pp) {
    pp <- pp +
        labs(x = expression(N['CUR']~(diploid)), y = expression(T['CUR']~(generation)))
    return(pp)
}

change_layer <- function(ppB) {
    ppB$layers[[5]] <- ppB$layers[[4]]
    ppB$layers[[4]] <- ppB$layers[[3]]
    ppB$layers[[3]] <- ppB$layers[[2]]
    ppB$layers[[2]] <- ppB$layers[[5]]
    ppB$layers[[5]] <- NULL
    return(ppB)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

outdir = './manuscript_plots/'
dir.create(outdir)

strokesize = 1

# load data --------
gridbig = readRDS('./ENP_1D.3Epoch_20210804/derived_data/dadi.grid.search.ENP.1D.3Epoch.LL.Rprocessed.rds')
griddiag = readRDS('./ENP_1D.3Epoch_2LL_20210804/derived_data/dadi.grid.search.ridge.ENP.1D.3Epoch.LL.Rprocessed.rds')
# x = nuF_scaled_dip_fixtheta, y = TF_scaled_gen_fixtheta
ppA = readRDS('./ENP_1D.3Epoch_20210804/plots/ENP.1D.3Epoch.gridsearch_scaled_fixtheta_20210805.rds')
ppB = readRDS('./ENP_1D.3Epoch_20210804/plots/ENP.1D.3Epoch.gridsearch_dLL_ann2_scaled_fixtheta_20210805.rds')
ppC = readRDS('./ENP_1D.3Epoch_2LL_20210804/plots/ENP.1D.3Epoch.gridridge_dLL_ann_scaled_fixtheta_20210805.rds')
# x = nuF_scaled_dip_fixtheta, y = TF_scaled_gen_fixtheta
ppD = readRDS('./ENP_1D.3Epoch_2LL_20210804/plots/ENP.1D.3Epoch.gridridge_dLL_ann2_scaled_fixtheta_20210805.rds')


# main --------
ppA$layers[[2]]$aes_params$size <- 3
ppA$layers[[3]]$aes_params$size <- 3
ppA$layers[[2]]$aes_params$stroke <- strokesize
ppA$layers[[3]]$aes_params$stroke <- strokesize

ppB$layers[[2]]$aes_params$stroke <- strokesize
ppB$layers[[3]]$aes_params$stroke <- strokesize
ppB <- change_layer(ppB)
ppB <- ppB +
    geom_mark_hull(data = griddiag, expand = 0, radius = 0, fill = NA)

ppD <- ppD +
    scale_fill_gradientn(breaks=c(-30,-2,-1,-0.1,-0.001),limits=c(-30,0),colors =c("darkgray","darkorange","yellow","green"),labels = c('-30','-2','-1','-0.1','-0.001'),trans = 'exp')
ppD$layers[[2]]$aes_params$stroke <- strokesize
ppD$layers[[3]]$aes_params$stroke <- strokesize
ppD$layers[[4]]$aes_params$size <- 3
ppD$layers[[4]]$aes_params$stroke <- 0.3
ppD <- change_layer(ppD)

pplist <- list(ppA, ppB, ppD)
pplist <- lapply(pplist, change_labs)
outpp <- ggarrange(plotlist = pplist, nrow = 3, ncol = 1, labels = c('A','B','C'))

ggsave(filename = paste0(outdir, 'ENP.1D.3Epoch.gridsearch_suppplot_', today, '.pdf'), plot = outpp, height = 12, width = 6)
# ggsave(filename = paste0(outdir, 'ENP.1D.3Epoch.gridsearch_suppplot_20210805.png'), plot = outpp, height = 8, width = 10)

# output files --------
# Modification: Add source_data
# Date: Mon Jan 23 10:49:13 2023
# ppA and ppB are from the same data sources
ppA_data = ppA$data %>%
    dplyr::select(nuB, nuF, TB, TF, nuF_scaled_dip_fixtheta, TF_scaled_gen_fixtheta,
                  LL_model, LL_data, LL_diff)

ppD_data = ppD$data %>%
    dplyr::select(nuB, nuF, TB, TF, nuF_scaled_dip_fixtheta, TF_scaled_gen_fixtheta,
                  LL_model, LL_data, LL_diff)
write.csv(ppA_data, file = '~/Lab/fin_whale/FinWhale_PopGenomics_2021/source_data/FigS14ab.csv')
write.csv(ppD_data, file = '~/Lab/fin_whale/FinWhale_PopGenomics_2021/source_data/FigS14c.csv')

# cleanup --------
date()
closeAllConnections()
