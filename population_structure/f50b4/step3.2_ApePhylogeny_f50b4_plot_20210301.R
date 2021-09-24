# Title: Plot the ape neighbor joining tree
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Mon Mar  1 02:15:28 2021
# Example: Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/f50b4/step3.2_ApePhylogeny_f50b4_plot_20210301.R '/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke' 'f50b4_pass_bialleic_all_LDPruned_maf05'

# IMPORTANT NOTE: The input rows (sample names) for ape::nj need to be named with values other than 1,2,3... The dist.gene function will confuse the dimmension of genotype
# IMPORTANT NOTE: Here the input matrix should be on the dosage of alternative alleles (when you use the dosage of reference allele, tree looks differently.)
# IMPORTANT NOTE: The ggtree is pretty outdated in bioconda. https://github.com/YuLab-SMU/ggtree/issues/91

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(gdsfmt)
library(SNPRelate)
# library(SeqArray)
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)

source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R')
# source('/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/plotting_config.R')

# def functions --------
# get maf
get_maf <- function(outprefix) {
    maf = stringr::str_split(outprefix, pattern = '_')[[1]]
    maf = maf[length(maf)]
    return(maf)
}

# empty legend
theme_blank_legend <- theme(
    legend.position = 'bottom',
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent", colour = NA), # get rid of legend panel bg
    legend.key = element_rect(fill = "transparent", colour = NA) # get rid of key legend fill, and of the surrounding
)

plot_trees <- function(tree_nj, type = c('rooted', 'unrooted'), popmap, outgrp, main, nodecutoff = 750) {
    metadt0 <- dplyr::left_join(data.frame('label' =  tree_nj$tip.label), popmap,
                               by = c('label' = 'SampleId'))
    if(type == 'unrooted') {
        p1_un <- ggtree(tree_nj)
        p1dt <- p1_un$data
        metadt <- dplyr::full_join(p1dt, metadt0, by = 'label') %>%
            dplyr::mutate(nodebt = ifelse(isTip == FALSE, label, NA))
        p1_un$data <- metadt
        outpp <- p1_un +
            geom_nodepoint(aes(size = as.integer(nodebt)), alpha = 0.5) +
            geom_tiplab(aes(color = SubPopId), show.legend = FALSE) +
            scale_color_manual(values = c(speccolors[-5],loccolors)) +
            labs(title = main, size = 'Bootstrap Support') +
            coord_cartesian(xlim = c(min(metadt$x), 1.1*max(metadt$x))) +
            theme_blank_legend
    } else {
        tree_nj_root <- ape::root(tree_nj, outgroup = outgrp, resolve.root = TRUE, interactive = FALSE, edgelabel = FALSE)
        p1_rt <- ggtree(tree_nj_root)
        p1dt <- p1_rt$data
        metadt <- dplyr::full_join(p1dt, metadt0, by = 'label') %>%
            dplyr::mutate(nodebt = ifelse(isTip == FALSE & label != 'Root', label, NA)) %>%
            dplyr::mutate(nodebt = ifelse(as.integer(nodebt) >= 750, nodebt, NA))
        p1_rt$data <- metadt
        outpp <- p1_rt +
            geom_tiplab(aes(color = SubPopId), show.legend = FALSE) +
            geom_point(aes(size = as.integer(nodebt), fill = as.integer(nodebt)), shape = 21, alpha = 0.8) +
            # geom_label(aes(label = nodebt), nudge_x = 1, size = 3, alpha=0.5) +
            scale_size_binned(range = c(1,8)) + 
            scale_fill_gradient(low = 'white', high = 'red') +
            scale_color_manual(values = c(speccolors[-5],loccolors)) +
            labs(title = main, size = 'Bootstrap Support', fill = '') +
            coord_cartesian(xlim = c(min(metadt$x), 1.1*max(metadt$x))) +
            theme_blank_legend
    }
    return(outpp)
}

# plot the density trees
plot_densitrees <- function(btrees, alpha = 0.1, main) {
    pp <- ggdensitree(btrees, alpha = alpha, colour='steelblue') +
        geom_tiplab(size=3) +
        labs(title = main)
    return(pp)
}

# def variables --------
args <- commandArgs(trailingOnly=TRUE)
workdir <- as.character(args[1]) # the working directory (should be the same as before and after)
outprefix <- as.character(args[2]) # the output prefix

today = format(Sys.Date(), "%Y%m%d")
datadate = '20210302'
outgrp = 'EubGla01'
nbt = 1000
maf = get_maf(outprefix)
outdir = './derive_phylogeny/'
plotdir = './plots/'

setwd(workdir)
sessionInfo()
# dir.create(outdir)
dir.create(plotdir)

# load data --------
# Usually, branch lengths can be interpreted as an estimate for the substitution
treespair <- readRDS(file = paste0(outdir, outprefix, '_treespair_20210302.rds'))
treesperc <- readRDS(file = paste0(outdir, outprefix, '_treesperc_20210302.rds'))

# get the population names
popmap = read.csv(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/baleen_popid.csv', stringsAsFactor = FALSE)
# popmap = read.csv(file = '/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/baleen_popid.csv', stringsAsFactor = FALSE)

# main --------

# plot the ape trees ========
pp_pair <- plot_trees(tree_nj = treespair[[1]], type = 'rooted', popmap, outgrp, main = paste0('NJ pairwise Root:EubGla ', maf))
pp_perc <- plot_trees(tree_nj = treesperc[[1]], type = 'rooted', popmap, outgrp, main = paste0('NJ percentage Root:EubGla ', maf))

ggsave(filename = paste0(outprefix, '_RtreePair_', today, '.pdf'), path = plotdir, plot = pp_pair, width = 8, height = 8)

ggsave(filename = paste0(outprefix, '_RtreePerc_', today, '.pdf'), path = plotdir, plot = pp_perc, width = 8, height = 8)

# # plot the density trees ========
# p2_pair <- plot_densitrees(btrees = treespair[[2]], main = paste0('Bootstrap trees pairwise Root:EubGla ', maf))
# p2_perc <- plot_densitrees(btrees = treesperc[[2]], main = paste0('Bootstrap trees percentage Root:EubGla ', maf))

# ggsave(filename = paste0(outprefix, '_BtreePair_', today, '.png'), path = plotdir, plot = p2_pair, width = 8, height = 8)

# ggsave(filename = paste0(outprefix, '_BtreePerc_', today, '.png'), path = plotdir, plot = p2_perc, width = 8, height = 8)

# cleanup --------
date()
closeAllConnections()
