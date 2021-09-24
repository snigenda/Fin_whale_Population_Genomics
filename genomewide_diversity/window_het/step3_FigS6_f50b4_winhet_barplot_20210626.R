# Title: FINAL Figure S6 Window-based Heterozygosity plots for baleen whale and GOC whale comparisons
# Author: Jacqueline Robinson, Paulina Nunez Valencia (pnunez@lcg.unam.mx), Sergio Nigenda and Meixi Lin (meixilin@ucla.edu)
# Date: Sat Sep 11 18:37:15 2021
# NOTE: Here I used the all50 dataset for the fin whales and pulled data from the f50b4 dataset
# Date: Mon Jun  7 18:14:42 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

source("~/Lab/finwhale_manuscript/scripts/config/plotting_config.R")

# def functions --------
# format genomewide het to annotate 
read_genomehet <- function(iucndt) {
    all50het = read.csv('./data/genome_stats/derived_data/all50_genomewide_heterozygosity_20210824.csv') %>%
        dplyr::filter(SampleId == 'GOC002')
    f50b4het = read.csv('./data/genome_stats/derived_data/f50b4_genomewide_heterozygosity_20210824.csv') %>%
        dplyr::filter(SampleId %in% c('BalAcu02', 'MegNov01', 'BalMus01', 'EubGla01'))
    outhet = left_join(iucndt, dplyr::bind_rows(all50het, f50b4het), by = 'SampleId')
    outhet = outhet %>% 
        dplyr::mutate(hetpkb = stringr::str_sub(as.character(1000*GenomeHet), end = 5)) %>% 
        dplyr::mutate(anntext = paste('Mean het. =', hetpkb, 'IUCN Status = ', IUCN)) %>%
        dplyr::select(SampleId, anntext) %>% 
        dplyr::rename(sample = SampleId)
    outhet$sample = factor(outhet$sample, levels = iucndt$SampleId)
    return(outhet)
}

# raw histogram method for any given vectors (left closed)
raw_hist <- function(hetpkb, breaks) {
    hetcat = vector(length = length(breaks)-1)
    for (ii in 1:length(breaks)-1) {
        jj = ii + 1
        hetcat[ii] = sum(hetpkb >= breaks[ii] & hetpkb < breaks[jj])
    }
    return(hetcat)
}

# nwin_all50 = 2284; nwin_f50b4 = 2258
get_histdt0 <- function(plotdt, histbreaks) {
    sample = unique(plotdt$sample)
    nsample = length(sample)
    nbins = length(histbreaks)-1
    histdt = data.frame(matrix(ncol = nsample+2, nrow = nbins))
    colnames(histdt) = c('minhetpkb_in', 'maxhetpkb_ex', sample)
    histdt$minhetpkb_in = histbreaks[-length(histbreaks)]
    histdt$maxhetpkb_ex = histbreaks[-1]
    for (ii in sample) {
        hetpkb = plotdt[plotdt$sample == ii,'hetpkb']
        hetcat = raw_hist(hetpkb = hetpkb, breaks = histbreaks)
        # summary of hetcat should be nwin
        histdt[,ii] = hetcat
    }
    return(histdt)
}

# tally up bottom of a histogram (note that all these are left closed)
tally_bottom <- function(histdt, cutoff) {
    dt1 = histdt[histdt$maxhetpkb_ex <= cutoff,] %>% 
        dplyr::mutate(maxhetpkb_ex = as.character(maxhetpkb_ex))
    dt2 = histdt[histdt$maxhetpkb_ex > cutoff,] 
    sumdt = data.frame(t(colSums(dt2)))
    sumdt$maxhetpkb_ex = paste0('>=', cutoff)
    outdt = rbind(dt1, sumdt) %>%
        dplyr::select(-minhetpkb_in) 
    # calculate the percents
    for (ii in 2:ncol(outdt)) {
        outdt[,ii] <- outdt[,ii]/sum(outdt[,ii])
    }
    print(colSums(outdt[,-1]))
    return(outdt)
}

# def variables --------
dataset = "f50b4"
nsample = 54
ref = "Minke"
today = format(Sys.Date(), "%Y%m%d")
datadate = "20210128" # the date the data generated

# window het settings
winsize = 1e+6
stepsize = 1e+6
misscutoff = 0.2 # minimum called/winsize (same as isle royale paper)

# IUCN list
SampleId = c('BalAcu02', 'MegNov01', 'BalMus01', 'EubGla01', 'GOC002')
IUCN = c('LC','LC','EN','EN','VU')
iucndt = data.frame(SampleId, IUCN)

# set axis limits
axis_lim = c(0,9.5,0,650)

# go to that dataset
indir = './data/window_het/derived_data/'
outdir = './data/window_het/derived_data/'
plotdir = './plots/window_het/'
dir.create(plotdir)
sessionInfo()

# load data --------
hetanndt <- read_genomehet(iucndt)
all50plotdt <- readRDS(file = './data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_20210824.rds')
f50b4plotdt <- readRDS(file = './data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_f50b4_plotdt_20210824.rds')

goc002 = all50plotdt %>% dplyr::filter(sample == 'GOC002')
baleen = f50b4plotdt %>% dplyr::filter(sample %in%  c('BalAcu02', 'MegNov01', 'BalMus01', 'EubGla01'))
plotdt1 = bind_rows(baleen, goc002)
# order by IUCN status
plotdt1$sample = factor(plotdt1$sample, levels = SampleId)

# main -------
# get the colors
mypalette = c(speccolors, loccolors)[unique(plotdt1$subpop)]
# get axis limit
max(plotdt1$hetpkb)
# window het ========
# the size option changes the width of the column
pp <- ggplot() +
    geom_col(data = plotdt1, aes(x = merge_start, y = hetpkb, color = subpop, fill = subpop), size = 1e-1) +
    geom_text(data = hetanndt, aes(x = 1e+8, y = 8, label = anntext), hjust = 'left') +
    scale_color_manual(values = mypalette) +
    scale_fill_manual(values = mypalette) +
    labs(y="Het./kb", x= "Scaffolds (bp)") +
    facet_wrap(. ~ sample, ncol = 1) +
    coord_cartesian(ylim=axis_lim[1:2], expand = FALSE) +
    ggpubr::theme_pubr() +
    theme(# axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none',
          strip.background = NULL)

# histograms ========
# changed the closed to left, but seemed fine both ways
pph <- ggplot(plotdt1, aes(x = hetpkb, color = subpop, fill = subpop)) +
    geom_histogram(binwidth = 0.15, size = 0.01, closed = 'left') +
    facet_wrap(. ~ sample, ncol = 1) +
    scale_color_manual(values = mypalette) +
    scale_fill_manual(values = mypalette) +
    labs(y="# of Windows", x= "Het./kb") +
    coord_cartesian(xlim=axis_lim[1:2], ylim=axis_lim[3:4], expand = TRUE) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none',
          strip.background = NULL)

# get a barplot that plots all the GOC individuals ========
# use the all50 GOC data (avoid confusion with figure 2)
gocdt = all50plotdt %>% dplyr::filter(subpop == 'GOC')
gocsamples = unique(gocdt$sample)
length(unique(gocdt$sample))

plotdt2 = rbind(baleen, gocdt)
unique(plotdt2$sample)
# order by IUCN status
if (length(c('BalAcu02', 'MegNov01', 'BalMus01', 'EubGla01', gocsamples)) != length(unique(plotdt2$sample))) {
    stop('Wrong plotdt2 label')
}

# generate data for barplot ========
# tested that previous ways of generating the histogram from ggplot was too complicated
# use a raw method (same as step4)
rawhistdt <- get_histdt0(plotdt2, histbreaks = seq(axis_lim[1], axis_lim[2], by = 0.1))
print(colSums(rawhistdt))
# tally up the bottom
barhistdt <- tally_bottom(histdt = rawhistdt, cutoff = 1)
barhistdt$maxhetpkb_ex = factor(barhistdt$maxhetpkb_ex, levels = rev(barhistdt$maxhetpkb_ex))
barplotdt <- reshape2::melt(barhistdt, id.var = 'maxhetpkb_ex', variable.name = 'sample', value.name = 'percent') 
barplotdt$sample = factor(barplotdt$sample, levels = c('BalAcu02', 'MegNov01', 'BalMus01', 'EubGla01', gocsamples))

# plot barplot of histogram ========
ppbar <- ggplot(data = barplotdt, aes(x = sample, y = percent, fill = maxhetpkb_ex)) +
    geom_bar(stat = 'identity') +
    scale_fill_viridis_d() + 
    scale_y_continuous(labels = scales::percent) + 
    labs(y = "% of 1Mb Windows", fill = "Het/kb") +
    ggpubr::theme_pubr() +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          legend.position = 'right')

# arrange them together ========
ppout <- ggarrange(ggarrange(pp, pph, widths = c(3,1), labels = c('A','B')),
                   ppbar,
                   nrow = 2, heights = c(2,1.2), labels = c('','C'))

ggsave(filename = paste0('FigureS6.winHet_1Mbwin_1Mbstep_20Per_three_f50b4_', today,'.pdf'), plot = ppout, path = plotdir,
       height = 12, width = 9)

# output data --------
write.csv(barhistdt, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_4baleen20goc_barhistdt_', today, '.csv'))

# cleanup --------
closeAllConnections()
date()

