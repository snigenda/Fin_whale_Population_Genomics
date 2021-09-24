# Title: FINAL Plot the Rxy for all50_snpEff_matching_PASSm6 dataset
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Apr 22 17:33:31 2021
# Modification: Use the new color schematics and remove the violin plot
# Date: Sun Sep  5 19:06:11 2021



# preparation --------
rm(list = ls())
cat("\014")

library(dplyr)
library(ggplot2)

# def functions --------

# plot the rxy
plot_rxy <- function(jk_forplot, rxy_summary, r2xy_summary, muttypes, mutlabs) {
    forplot = jk_forplot %>%
        dplyr::filter(runNum > 0 & MutType %in% muttypes) %>%
        dplyr::select(MutType, RRxy, RR2xy) %>%
        reshape2::melt(id.vars = 'MutType')
    summary2.1 = rxy_summary %>%
        dplyr::filter(MutType %in% muttypes) %>%
        dplyr::mutate(variable = 'Rxy')
    summary2.2 = r2xy_summary %>%
        dplyr::filter(MutType %in% muttypes) %>%
        dplyr::mutate(variable = 'R2xy')
    colnames(summary2.2) = colnames(summary2.1)

    summary = dplyr::bind_rows(summary2.1, summary2.2) %>%
        dplyr::mutate(z_plab = case_when(z_pvalue < 0.001 ~ '***',
                                         z_pvalue < 0.01 & z_pvalue >= 0.001 ~ '**',
                                         z_pvalue < 0.05 & z_pvalue >= 0.01 ~ '*',
                                         TRUE ~ 'ns'))
    summary$MutType = factor(summary$MutType, levels = rev(muttypes), labels = rev(mutlabs))

    forplot$MutType = factor(forplot$MutType, levels = rev(muttypes), labels = rev(mutlabs))
    forplot$variable = factor(forplot$variable, levels = c('RRxy', 'RR2xy'), labels = c('Rxy', 'R2xy'))
    pp <- ggplot() +
        geom_hline(yintercept = 1, color = 'darkgray', linetype = 'dashed') +
        # geom_violin(data = forplot, aes(x = MutType, y = value, color = MutType)) +
        geom_errorbar(data = summary, aes(x = MutType, ymin = low2sig_Rxy, ymax = up2sig_Rxy, color = MutType), width = 0.2) +
        geom_point(data = summary, aes(x = MutType, y = Rxy, color = MutType), size = 3) +
        geom_text(data = summary, aes(x = MutType, y = up2sig_Rxy+0.05, label = z_plab)) +
        facet_wrap(. ~  variable) +
        scale_color_manual(values = mycolors) +
        coord_flip() +
        ggpubr::theme_pubr() +
        labs(color = '', x = '', y = 'GOC/ENP Ratio')
    return(pp)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50_snpEff_matching'
ref = 'Minke'
gttype = 'PASSm6'
prefixlist = c('LOW','MODERATE','HIGH')

xpop = 'GOC'
ypop = 'ENP'

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/Rxy', dataset, ref, sep = '/')
setwd(workdir)

plotdir = './plots/'

sessionInfo()

mycolors = c('#C94571', '#838936', '#CD960C')
names(mycolors) = prefixlist

# load data --------
rxy_summary = read.csv(file = 'derive_data/rxy_table/Summary_all50_snpEff_matching_RxyP_xGOCyENP_PASSm6_20210422.csv',
                       stringsAsFactors = FALSE, row.names = 1)
r2xy_summary = read.csv(file = 'derive_data/rxy_table/Summary_all50_snpEff_matching_R2xyP_xGOCyENP_PASSm6_20210422.csv',
                        stringsAsFactors = FALSE, row.names = 1)
jk_forplot = readRDS("derive_data/rxy_table/Summary_all50_snpEff_matching_Minke_1000JK_xGOCyENP_PASSm6_20210422.rds")

# main --------

# plotting the jackknife and se estimates ========
pp1 <- plot_rxy(jk_forplot, rxy_summary, r2xy_summary, prefixlist, prefixlist)
ggsave(filename = paste0('NewColRxyImp_all50_snpEff_matching_Minke_1000JK_xGOCyENP_PASSm6_', today,'.pdf'),
       plot = pp1, path = plotdir, height = 4, width = 7)

# cleanup --------
date()
closeAllConnections()
