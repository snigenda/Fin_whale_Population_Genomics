# Title: Plot admixture
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Jan 11 17:33:58 2021
# Modification: Update to the new LDpruning schemes
# Date: Fri Mar 19 00:20:21 2021
# Modification: Generate source data
# Date: Mon Jan  9 13:32:32 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/Users/linmeixi/google_drive/finwhale/analyses/important_results/Manuscript_plots")
today = format(Sys.Date(), "%Y%m%d")

# sink(file = paste0("logs/admixture_cv_", today,".log"))
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# def functions --------

# def variables --------
mafcut='10'
# mafcut='05'
# mafcut='NA'

sessionInfo()
source("/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# load data --------
# read cross validation files
adcv = read.csv(file = paste0("/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/PopStructure/all50/Minke/Admixture_20210318/maf", mafcut, "/Admixture_CV_LLsummary_maf", mafcut, ".csv"), stringsAsFactors = FALSE)

forplot = adcv %>%
    dplyr::group_by(K) %>%
    dplyr::summarise(meancv = mean(CVERROR),
                     sdcv = sd(CVERROR))

# main --------
pp1 <- ggplot(forplot, aes(x = K, y = meancv)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=meancv-sdcv, ymax=meancv+sdcv), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = "K", y = "Cross-Validation Error") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 11))

ggsave(filename = paste0("admixture_maf", mafcut, "_cv_", today, ".pdf"), plot = pp1, height = 3, width = 3, bg = "transparent")

write.csv(forplot, file = '~/Lab/fin_whale/FinWhale_PopGenomics_2021/source_data/FigS3b.csv')

# cleanup --------
# sink()
closeAllConnections()
