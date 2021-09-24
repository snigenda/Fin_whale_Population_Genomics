# Title: FINAL manuscript plots for Fig 5 slim simulations
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Sep 10 23:30:05 2021
# Modification: Use the most uptodate ENP models
# Date: Sun Sep 19 22:55:50 2021


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
        # stat_summary(fun = "mean", geom = 'point', shape = 23, mapping = aes(fill = pop)) +
        scale_fill_manual(values = mycolors) +
        scale_color_manual(values = mycolors) +
        ggpubr::theme_pubr() +
        labs(y = ylabel, x = "") +
        theme(legend.position = 'none',
              text = element_text(family = 'ArialMT',
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
plotdir = './plots/maintext/'
indir = './data/slim/derived_data/'

# variables to collect
plotting_vars = c('meanHet', 'FROH_1Mb', 'geneticLoad')
ylabs = c('Heterozygosity', 'Inbreeding coefficient', 'Genetic load')

# load data --------
plot2popdt = readRDS(file = paste0(indir, 'slim_plotdt_2pop_ancestralChange_20210919.rds'))
plotenpdt = readRDS(file = paste0(indir, 'slim_plotdt_ENP_20210919.rds'))

# main --------
# panel A: enp ========
longenpdt = plotenpdt %>%
    reshape2::melt(id.vars = c('rep','model', 'pop')) %>% 
    dplyr::filter(variable %in% plotting_vars) %>%
    dplyr::mutate(popmodel = model)
longenpdt$variable = factor(longenpdt$variable, levels = plotting_vars)

# plotting ========
ppAlist <- vector('list',length = length(plotting_vars))
ppAlist[[1]] <- plot_set(longdt = longenpdt, plotvar = plotting_vars[1], 
                         ylabel = ylabs[1], ybreaks = seq(from = 0.0013, to = 0.0015, by = 0.0001))
ppAlist[[2]] <- plot_set(longdt = longenpdt, plotvar = plotting_vars[2], ylabel = ylabs[2])
ppAlist[[3]] <- plot_set(longdt = longenpdt, plotvar = plotting_vars[3], ylabel = ylabs[3]) 

ppAlist <- lapply(ppAlist, function(pp) {
    pp <- pp +
        scale_x_discrete(labels = c('pre-bott','2 gens','20 gens','20 gens\nw/ recov'))
})

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
outpp = ggpubr::ggarrange(ppAlist[[1]], ppBlist[[1]],
                          ppAlist[[2]], ppBlist[[2]],
                          ppAlist[[3]], ppBlist[[3]],
                          nrow = 3, ncol = 2, align = 'hv')

# output files --------
ggsave(filename = paste0(plotdir, 'Fig5_slimv2.pdf'), plot = outpp, width = 8.5, height = 9)

# cleanup --------
date()
closeAllConnections()

