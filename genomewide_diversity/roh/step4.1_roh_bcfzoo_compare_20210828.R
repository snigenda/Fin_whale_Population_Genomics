# Title: Output summaries for bcftools and zooroh comparisons
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Mar 25 16:24:29 2021
# Modification: Update to match the new categoriess
# Date: Thu Aug 12 16:55:25 2021
# Modification: Final supplemental plots (with name of categories)
# Date: Sat Aug 28 11:36:14 2021


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/Users/linmeixi/google_drive/finwhale/analyses/important_results/Runs_of_homozygosity/")

library(dplyr)
library(ggplot2)

# def functions --------
# generate rohcategory names
# rohlens should be in the units of Mb
get_rohnames <- function(rohlens) {
    output=character(0)
    for (ii in 1:length(rohlens)) {
        if (ii == length(rohlens)) {
            temp = paste0(rohlens[ii], '_Inf')
        } else {
            temp = paste0(rohlens[ii], '_',rohlens[ii+1])
        }
        output = c(output, temp)
    }
    return(output)
}

# get bcf and zooroh compare plot
get_forplot <- function(rohsum) {
    forplot = rohsum %>%
        select(-starts_with('froh'), -GenomeHet) %>%  
        reshape2::melt(., id.vars = c('SampleId', 'PopId', 'SubPopId')) 
    dt = reshape2::colsplit(forplot$variable, pattern = '_', names = c('type', 'software', 'length'))
    forplot = cbind(forplot, dt) 
    forplotbcf = forplot %>% 
        dplyr::filter(software == 'bcf') %>% 
        dplyr::rename(bcftools = value) %>%
        dplyr::select(SampleId, PopId, SubPopId, type, length, bcftools)
    forplotzoo = forplot %>% 
        dplyr::filter(software == 'zoo') %>% 
        dplyr::rename(RZooRoH = value) %>%
        dplyr::select(SampleId, PopId, SubPopId, type, length, RZooRoH)
    forplot = dplyr::full_join(x = forplotbcf, y = forplotzoo, 
                               by = c('SampleId','PopId','SubPopId','type','length')) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(type = ifelse(type == 'N', 'Total number', 'Total length (Mb)'))
    return(forplot)
}

# def variables --------
plotdir = "./plots/"
outdir = "./derived_data/"
dir.create(plotdir, recursive = T)
dir.create(outdir, recursive = T)

today = format(Sys.Date(), "%Y%m%d")
genomelen = 2324429847

# new roh length and categories (August 2021)
rohlens = c(0.1,1,5)*1e+6
rohcats = get_rohnames(c(0.1,1,5))

sessionInfo()
source("/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# load data --------
# load the newly calculated roh categories
rohsum2 = readRDS(file = 'derived_data/rohsummary_zoobcf_20210808.rds')

# plot comparisons ========
forplot <- get_forplot(rohsum2)
types <- unique(forplot$type)
lens <- unique(forplot$length)
# get labels for different lens
lenslab <- c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb', '[0.1, Inf) Mb')
lenscat <- c('short', 'medium', 'long', 'all')
plotlist <- vector(mode = "list", length = length(types)*length(lens))
counter = 0
for (ii in types) {
    for (j in 1:length(lens)) {
        jj = lens[j]
        lenlab = lenslab[j]
        lencat = lenscat[j]
        counter = counter + 1
        forplotx = forplot %>% 
            dplyr::filter(type == ii, length == jj) 
        if (ii == "Total length (Mb)") {
            forplotx = forplotx %>%
                dplyr::mutate(bcftools = bcftools/1e+6,
                              RZooRoH = RZooRoH/1e+6)
        }
        pltrange = range(forplotx$bcftools,forplotx$RZooRoH)
        my.formula = y ~ x
        pp <- ggplot() + 
            geom_point(data = forplotx, mapping = aes(x = bcftools, y = RZooRoH), shape = '.') + 
            geom_smooth(data = forplotx, mapping = aes(x = bcftools, y = RZooRoH), method = "lm", se=FALSE, formula = my.formula, size = 0.5, linetype = 'dashed', color = 'black') +
            ggpmisc::stat_poly_eq(data = forplotx, 
                                  mapping = aes(x = bcftools, y = RZooRoH, 
                                                label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                                  formula = my.formula, parse = TRUE, size = 3) +
            geom_smooth(data = forplotx, mapping = aes(x = bcftools, y = RZooRoH, color = PopId), 
                        method = "lm", se=FALSE, formula = my.formula, size = 0.8) +
            geom_point(data = forplotx, mapping = aes(x = bcftools, y = RZooRoH, color = PopId)) +
            ggpmisc::stat_poly_eq(data = forplotx, 
                                  mapping = aes(x = bcftools, y = RZooRoH,  color = PopId,
                                                label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                                  formula = my.formula, parse = TRUE, size = 3, label.y = c(0.9, 0.85)) +
            geom_abline(slope = 1, intercept = 0, color = 'darkgray', linetype = 'dotted') +
            labs(title = paste(ii, 'of', lencat, 'ROH'), subtitle = paste('ROH length:', lenlab)) + 
            coord_fixed(ratio = 1, xlim = pltrange, ylim = pltrange) +
            scale_color_manual(values = mycolors) +
            theme_bw() +
            theme(legend.position = 'none',
                  plot.title = element_text(size = 12),
                  plot.subtitle = element_text(size = 10))
        # The overlapping one was ENPCA01, expected
        plotlist[[counter]] <- pp
    }
}

ppout <- ggpubr::ggarrange(plotlist = plotlist, nrow = length(types), ncol = length(lens))

# output --------
ggsave(filename = paste0('FigureS8.ROH_bcfzoo_compare_', today, '.pdf'), path = plotdir, width = 14, height = 8)

# cleanup --------
date()
closeAllConnections()
