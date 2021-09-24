# Title: Plot the Rxy
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Mar  7 17:50:56 2021
# Modification: Update the new equation for R2xy
# Date: Fri Mar 19 18:15:31 2021
# Modification: Add new color scale
# Date: Fri Aug 13 18:11:36 2021
# Modification: Change font
# Date: Fri Sep  3 15:36:39 2021



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
        geom_errorbar(data = summary, aes(x = MutType, ymin = low2sig_Rxy, ymax = up2sig_Rxy, color = MutType), width = 0.2) +
        # geom_violin(data = forplot, aes(x = MutType, y = value, color = MutType)) +
        geom_point(data = summary, aes(x = MutType, y = Rxy, color = MutType), size = 3) +
        geom_text(data = summary, aes(x = MutType, y = up2sig_Rxy+0.08, label = z_plab)) +
        facet_wrap(. ~  variable) +
        scale_color_manual(values = mycolors) +
        coord_flip() +
        ggpubr::theme_pubr() +
        labs(color = '', x = '', y = 'GOC/ENP Ratio')
    return(pp)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
cdstype = 'ALLregions'
prefixlist1 = c("syn", "nonsynTOL", "nonsynDEL", "LOF")
mutlab1 = c('SYN', 'TOL', 'DEL', 'LOF')

xpop = 'GOC'
ypop = 'ENP'

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/Rxy', dataset, ref, sep = '/')
setwd(workdir)

plotdir = './plots/'
outdir1 = './derive_data/lxy_table/'
outdir2 = './derive_data/rxy_jktable/'
outdir3 = './derive_data/rxy_table/'

sessionInfo()

mycolors = c('#9B696E', '#C94571', '#838936', '#CD960C')
names(mycolors) = mutlab1

# load data --------
# jackknife for plotting
jk_forplot = readRDS(paste0(outdir3, "Summary_all50_Minke_ALLregions_1000JK_xGOCyENP_PASSm6_20210319.rds"))

# the updated rxy_summary
rxy_summary = read.csv(file = './derive_data/rxy_table/Summary_all50_Minke_ALLregions_RxyP_xGOCyENP_PASSm6_20210319.csv')
r2xy_summary = read.csv(file = './derive_data/rxy_table/Summary_all50_Minke_ALLregions_R2xyP_xGOCyENP_PASSm6_20210319.csv')

# main --------
# plotting the jackknife and se estimates ========
pp1 <- plot_rxy(jk_forplot, rxy_summary, r2xy_summary, prefixlist1, mutlab1) +
    theme(text = element_text(family = 'ArialMT'))

ggsave(filename = paste0('NewColRxyMut_all50_Minke_ALLregions_1000JK_xGOCyENP_PASSm6_', today,'.pdf'),
       plot = pp1, path = plotdir, height = 4, width = 6)

# output files --------

# cleanup --------
date()
closeAllConnections()
