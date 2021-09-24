# Title: Window-based Heterozygosity plot objects
# Author: Jacqueline Robinson, Paulina Nunez Valencia (pnunez@lcg.unam.mx), Sergio Nigenda and Meixi Lin (meixilin@ucla.edu)
# Date: Mon Jun 21 15:04:31 2021

# Output: Figure S4. Genomewide distribution of heterozygosity in all individuals
# Output: Figure S5. Histograms of genomewide distribution of heterozygosity in all individuals
# Execution: Rscript --vanilla step3_FigS4S5_all50_window_het_20210824.R

# preparation --------
rm(list = ls())
cat("\014")

options(echo = TRUE)

library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)

source("~/Lab/finwhale_manuscript/scripts/config/plotting_config.R")

# def functions --------

# def variables --------
dataset = "all50"
nsample = 50
ref = "Minke"
today = format(Sys.Date(), "%Y%m%d")
datadate = "20201206" # the date the data generated

# window het settings
winsize = 1e+6
stepsize = 1e+6
misscutoff = 0.2 # minimum called/winsize (same as isle royale paper)

# go to that dataset
setwd('~/Lab/finwhale_manuscript/')
indir = './data/window_het/derived_data/'
outdir = './data/window_het/derived_data/'
plotdir = './plots/window_het'
dir.create(plotdir)
sessionInfo()

# load data --------
plotdt = readRDS(file = './data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_20210824.rds')
# set the axis limits from previous knowledge (different from the one used in the main text figure 2)
max(plotdt$hetpkb) # 8.143495 (used to set the axis_limit to include all the data)
axis_lim = c(0,8.2,0,500)

# load the genomewide diversity data as well to annotate
genomehet = read.csv(file = './data/genome_stats/derived_data/all50_genomewide_heterozygosity_20210824.csv', stringsAsFactors = FALSE, row.names = 1) %>%
    dplyr::select(SampleId, GenomeHet) %>% 
    dplyr::mutate(GenomeHetText = paste('Mean het. =', sprintf("%.3f", round(1000*GenomeHet, digits = 3)))) %>%
    dplyr::rename(sample = SampleId)

# main -------
# Figure S4. Genomewide distribution of heterozygosity in all individuals ========
# the size option changes the width of the column
pp <- ggplot() +
    geom_col(data = plotdt, aes(x = merge_start, y = hetpkb, color = subpop, fill = subpop), size = 1e-7) +
    geom_text(data = genomehet, aes(label = GenomeHetText), color = 'black', x = 1e+8, y = 8, hjust = 0, vjust = 1) +
    scale_color_manual(values = loccolors) +
    scale_fill_manual(values = loccolors) +
    labs(y="Het./kb", x= "Scaffolds (bp)") +
    facet_wrap(. ~ sample, nrow = 10, ncol = 5) +
    coord_cartesian(ylim=axis_lim[1:2], expand = FALSE) +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
ggsave(filename = paste0('FigureS4.winHet_1Mbwin_1Mbstep_20Per_all50_', today,'.pdf'), plot = pp, path = plotdir,
       height = 12, width = 15)

# Figure S5. Histograms of genomewide distribution of heterozygosity in all individuals ========
# Not getting the histogram data here because the binwidth set the minimum at -0.075
pph <- ggplot(plotdt, aes(x = hetpkb, color = subpop, fill = subpop)) +
    geom_histogram(binwidth = 0.15, size = 0.01) +
    facet_wrap(. ~ sample, nrow = 10, ncol = 5) +
    scale_color_manual(values = loccolors) +
    scale_fill_manual(values = loccolors) +
    labs(y="# of Windows", x= "Het./kb") +
    coord_cartesian(xlim=axis_lim[1:2], ylim=axis_lim[3:4], expand = TRUE) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none')
ggsave(filename = paste0('FigureS5.winHet_1Mbwin_1Mbstep_20Per_all50_hist_', today,'.pdf'), plot = pph, path = plotdir,
       height = 12, width = 9)

# output data --------

# cleanup --------
closeAllConnections()
date()
