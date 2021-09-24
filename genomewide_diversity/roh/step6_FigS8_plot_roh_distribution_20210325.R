# Title: Plot ZooRoH output along scaffolds
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Mar 25 20:26:32 2021
# Modification: Migrate to a new folder
# Date: Sun Sep 12 15:14:42 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

library(dplyr)
library(ggplot2)
library(ggpubr)

# def functions --------
load_contiglist <- function(contigfile) {
    contiglist = read.csv(file = contigfile, row.names = 1, stringsAsFactors = F)
    # change to zero based and only select the two columns relavant
    contiglist = contiglist %>%
        dplyr::mutate(gnPOS = genomewide_coord -1) %>%
        dplyr::select(SN, gnPOS)
    colnames(contiglist) = c('chrom', 'gnPOS')
    return(contiglist)
}

# def variables --------
# source("./scripts/config/plotting_config.R")
indir = "./data/roh/derived_data/"
outdir = "./data/roh/derived_data/"
plotdir = './plots/roh/'
dir.create(plotdir)

today = format(Sys.Date(), "%Y%m%d")
genomelen = 2324429847

# new roh length and categories (August 2021)
rohlens = c(0.1,1,5)*1e+6
rohcats = c('0.1_1','1_5','5_Inf')
lenslab <- c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb')

sessionInfo()

# load data --------
zooroh = readRDS(file = './data/roh/derived_data/ROH_zooroh_summary_final_20210907.rds')
contiglist = load_contiglist(contigfile = "./scripts/config/minke_contig_summary.csv")

# main --------
# plotting zooroh ========
nrohcat = length(unique(zooroh$rohcat)) # nrohcat =3 
rohcatbrewer = RColorBrewer::brewer.pal(name = "Reds", n = 9)[c(5,7,9)]
names(rohcatbrewer) = unique(zooroh$rohcat)
# pals::pal.bands(rohcatbrewer)
# https://stackoverflow.com/questions/28391850/reverse-order-of-discrete-y-axis-in-ggplot2
pp2 <- ggplot(zooroh) +
    geom_segment(aes(x=gnstart, xend=gnend,y=sample,yend=sample, color= rohcat), size = 3) +
    scale_color_manual(values = rohcatbrewer, labels = lenslab) +
    scale_y_discrete(limits = rev, expand = c(0,1)) +
    scale_x_continuous(limits = c(0, genomelen), expand = c(0,1)) +
    theme_bw() +
    labs(x = "ROH position in reference", y = "Individual", color = 'Category') +
    theme(legend.position = 'top')

# plotting only the long ones ========
plotzoo2 = zooroh %>% dplyr::filter(rohcat == '5_Inf')

pp3 <- ggplot(plotzoo2) +
    geom_vline(data = contiglist, aes(xintercept = gnPOS), color = 'darkgray', linetype = 'dotted')  +
    geom_segment(aes(x=gnstart, xend=gnend,y=sample,yend=sample, color= rohcat), size = 3) +
    scale_color_manual(values = rohcatbrewer) +
    scale_y_discrete(limits = rev, expand = c(0,1)) +
    scale_x_continuous(limits = c(0, genomelen), expand = c(0,1)) +
    theme_bw() +
    labs(x = "ROH position in reference", y = "Individual", color = 'Category') +
    theme(legend.position = 'none')

# arrange the two ========
pp <- ggarrange(pp2, pp3, nrow = 2, labels = 'AUTO', heights = c(2, 1), common.legend = TRUE)

ggsave(filename = 'FigureS8.genomewide_distribution_RZooRoH_newcat_20210912.pdf', plot = pp, path = plotdir, width = 12, height = 9)

# cleanup --------
closeAllConnections()