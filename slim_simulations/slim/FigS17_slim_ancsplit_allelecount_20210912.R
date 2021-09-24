# Title: FINAL manuscript plots for Fig S17 slim simulations allele count
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Sep 10 23:30:05 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

library(dplyr)
library(ggplot2)
library(ggpubr)

# def functions --------
plot_set <- function(longdt, plotvar, ylabel, ybreaks = integer(0)) {
    thisdt = longdt %>% filter(variable == plotvar)
    pp <- ggplot(data = thisdt, aes(x = popmodel, y = value)) +
        geom_boxplot(aes(color = pop, fill = pop), alpha = 0.5) +
        geom_point(color = 'gray10', shape = 0, stroke = 0.3) +
        # ggpubr::stat_compare_means(aes(label = ..p.signif..),
        #                            method = 'wilcox.test') +
        # stat_summary(fun = "mean", geom = 'point', shape = 23, mapping = aes(fill = pop)) +
        scale_fill_manual(values = mycolors) +
        scale_color_manual(values = mycolors) +
        ggpubr::theme_pubr() +
        labs(y = ylabel, x = "") +
        theme(legend.position = 'none',
              text = element_text(# family = 'ArialMT',
                                  size = 11))
    if (length(ybreaks) > 0) {
        pp <- pp +
            scale_y_continuous(breaks = ybreaks)
    }
    return(pp)
}

# def variables --------
source('~/Lab/fin_whale/scripts_analyses/config/plotting_config.R')
today = format(Sys.Date(), "%Y%m%d")
plotdir = './plots/slim/'
dir.create(plotdir)
indir = './data/slim/derived_data/'

# variables to collect
plotting_vars = c("avgStrDel","avgModDel","avgWkDel")
ylabs = c('# str.del. alleles', '# mod.del. alleles', '# weak del. alleles')

# load data --------
plot2popdt = readRDS(file = paste0(indir, 'slim_plotdt_2pop_ancestralChange_20210910.rds'))

# main --------
# panel B: 2 pop ========
long2popdt = plot2popdt %>%
    reshape2::melt(id.vars = c('rep','model', 'pop', 'popmodel')) %>% 
    dplyr::filter(variable %in% plotting_vars)
long2popdt$variable = factor(long2popdt$variable, levels = plotting_vars)

# plotting ========
ppBlist <- vector('list',length = length(plotting_vars))
ppBlist[[1]] <- plot_set(longdt = long2popdt, plotvar = plotting_vars[1], ylabel = ylabs[1])
ppBlist[[2]] <- plot_set(longdt = long2popdt, plotvar = plotting_vars[2], ylabel = ylabs[2])
ppBlist[[3]] <- plot_set(longdt = long2popdt, plotvar = plotting_vars[3], ylabel = ylabs[3])

# arrange =========
outpp = ggpubr::ggarrange(plotlist = ppBlist, nrow = 3, ncol = 1, align = 'hv')

# output files --------
ggsave(filename = 'FigureS17.Slim_del_var_2popAncChange_20210912.pdf', path = plotdir, plot = outpp, width = 5, height = 8)

# cleanup --------
date()
closeAllConnections()

