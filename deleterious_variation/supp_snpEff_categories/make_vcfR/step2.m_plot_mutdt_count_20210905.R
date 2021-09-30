# Title: Plot the master dataset with the invariant sites and bad individuals filtered out
# Note: SnpEff impact HIGH/MODERATE/LOW matching version
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Apr 22 16:37:56 2021
# Modification: FINAL use average counts
# Date: Sun Sep  5 18:26:48 2021


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

source("~/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# def functions --------
# from: step3.2_plot_mutdt_20210521.R
source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/facet_functions.R')
# from: step4_generate_summarystats_20210218.R in all50
# get mean normalized variations (by mean value)
get_mean_rawdt <- function(rawdt) {
    # from nopopdt, get Average CalledCount and CalledAlleleCount
    meandt = rawdt %>%
        dplyr::select(CalledCount, MutType) %>%
        group_by(MutType) %>%
        summarise(meanCalledCount = mean(CalledCount),
                  meanCalledAllele = 2*meanCalledCount,
                  .groups = 'drop')
    outdt = dplyr::left_join(x=rawdt, y = meandt, by = 'MutType') %>%
        dplyr::mutate(normHomRef = HomRefPer * meanCalledCount,
                      normHomAlt = HomAltPer * meanCalledCount,
                      normHet = HetPer * meanCalledCount,
                      normRefAllele = RefAllelePer * meanCalledAllele,
                      normAltAllele = AltAllelePer * meanCalledAllele)
    return(outdt)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50_snpEff_matching'
ref = 'Minke'
gttype = 'PASSm6'
# gttype = 'PASSm6CL' # IMPORTANT: DON'T USE THIS VERSION!!!
prefixlist = c('HIGH','MODERATE','LOW')

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR/', dataset, ref, sep = '/')
setwd(workdir)

plotdir = './pub_plots/'
dir.create(plotdir)

sink(file = paste0('logs/step3_plot_mutdt_count_', paste(dataset, ref, gttype, today, sep = '_'), ".log"))
date()
sessionInfo()

# load data --------
rawdt = readRDS(file = './derive_data/sum_table/HML_all50_snpEff_matching_Minke_SUMtable_Seg_PASSm6_20210422.rds')

# convert percentages to number of average variations ========
normdt = get_mean_rawdt(rawdt)

# output this
saveRDS(normdt, file = './derive_data/sum_table/HML_all50_snpEff_matching_Minke_SUMtable_Seg_PASSm6_Norm_20210905.rds')

# main --------
mutdt = normdt

mutdt$MutType = factor(mutdt$MutType, levels = prefixlist)

# the allele proportion plot ========
alleledt = reshape2::melt(mutdt, id.vars = c("SampleId", "PopId", "SubPopId", "MutType")) %>%
    dplyr::filter(variable %in% c("normAltAllele"))

pp1 <- ggplot(data = alleledt, aes(x = PopId, y = value, fill = PopId)) +
    geom_boxplot() + 
    scale_fill_manual(values = mycolors) +
    # facet_wrap_custom(. ~ MutType, scales = "free_y",scale_overrides = list(
    #     scale_override(1, scale_y_continuous(breaks = c(25400, 25600, 25800))),
    #     scale_override(2, scale_y_continuous(breaks = c(13700, 13900))),
    #     scale_override(3, scale_y_continuous(breaks = c(5800, 6000))),
    #     scale_override(4, scale_y_continuous(breaks = c(360, 380)))
    # )) +
    facet_wrap(. ~ MutType, scales = "free_y") +
    labs(fill = "Population", x = "Population", y = "Derived Allele Count") +
    theme_pubr() +
    theme(legend.position = "none")

pp1.2 <- pp1 +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox.test",
                       label.x.npc = 'center',
                       vjust = 1,
                       size = 4) 

ggsave(filename = paste0(plotdir, "SnpEffType_allele_number_box_", gttype, "_", today, ".pdf"), plot =  pp1.2, height = 4, width = 7)

# the genotype proportion plot ========
genodt = reshape2::melt(mutdt, id.vars = c("SampleId", "PopId", "SubPopId", "MutType")) %>% 
    dplyr::filter(variable %in% c("normHet", "normHomAlt"))
genodt$variable = as.character(genodt$variable)
genodt$variable = factor(genodt$variable, levels = c("normHomAlt", "normHet"), labels = c("Homozygous Derived", "Heterozygous"))

pp2 <- ggplot(data = genodt,
              aes(x = PopId, y = value, fill = PopId)) +
    geom_boxplot() + 
    scale_fill_manual(values = mycolors) +
    facet_grid(MutType ~ variable, scale = 'free_y') +
    # facet_grid_custom(MutType ~ variable, scale = 'free_y',scale_overrides = list(
    #     # scale_override(1, scale_y_continuous(limits = c(5500, 11500))),
    #     scale_override(2, scale_y_continuous(limits = c(3000, 6500))),
    #     scale_override(3, scale_y_continuous(limits = c(1200, 3200)))
    #     # scale_override(4, scale_y_continuous(breaks = c(100, 150)))
    # )) +
    labs(fill = "Population", x = "Population", y = "Genotype Count") +
    theme_pubr() +
    theme(legend.position = "none")

pp2.2 <- pp2 +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox.test",
                       label.x.npc = 'center',
                       vjust = 1,
                       size = 4) 


ggsave(filename = paste0(plotdir, "SnpEffType_geno_number_box_", gttype, "_", today, ".pdf"), plot =  pp2.2, height = 7, width = 5)

# cleanup --------
sink()
closeAllConnections()
