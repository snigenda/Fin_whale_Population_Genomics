# Title: Plot admixture output for all scenarios
# USED IN: SUPPLEMENTAL FIGURE 2
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Mar 22 23:56:00 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(ggplot2)
library(ggpubr)
library(dplyr)

source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R')

# def functions --------
read_admixture <- function(k, prefix, nrep, popmap, subpoporder) {
    # Q (the ancestry fractions) <-- What we need
    # P (the allele frequencies of the inferred ancestral populations).
    qfile = paste0(prefix, '.K', k,'.iter',1:nrep, '.Q')
    admaster = data.frame()
    for (ii in 1:length(qfile)) {
        ad = read.table(file = qfile[ii], stringsAsFactors = FALSE, header = FALSE)
        colnames(ad) = paste0("Cluster", 1:k) 
        ad$run = ii
        ad = base::cbind(popmap, ad)
        admaster = base::rbind(admaster, ad)
    }
    # append info on K
    admaster = admaster %>% dplyr::mutate(K = k)
    return(admaster)
}

format_admaster <- function(admaster, keeprun = c(2,4,6,8,10)) {
    K = unique(admaster$K)
    forplot = admaster %>%
        dplyr::filter(run %in% keeprun) %>%
        reshape2::melt(id.vars = c("SampleId", "PopId", "SubPopId", "K", "run"), measure.vars = paste0('Cluster', 1:K)) %>%
        dplyr::mutate(K = paste('K =', K))
    return(forplot)
}

get_sampleorder <- function(popmap, subpoporder) {
    popmap$SubPopId = factor(popmap$SubPopId, levels = subpoporder)
    dt = popmap %>%
        dplyr::arrange(SubPopId)
    sampleorder = unique(dt$SampleId)
    return(sampleorder)
}
# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
mafcut = '10'
Klist = 2:6
subpoporder = c("AK", "BC", "WA", "OR", "CA", "GOC") # order for subpopulations 

gdsfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds'
# confirm the sample names
plinkfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10_SA_mrF.nosex'

# set derived data outside of the 'Minke' folder to make hoffman2 folder unique
workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/PopStructure', dataset, sep = '/')
indir = paste0('./', ref, '/Admixture_20210318/')
outdir = './derive_data/'
plotdir = './plots/'

setwd(workdir)
sessionInfo()

dir.create(outdir)
dir.create(plotdir)

# load data --------
popmap = read.csv(file = "/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/popmap_all50.csv", stringsAsFactors = F)
plinkname = read.table(file = paste0('./Minke/', plinkfile), header = FALSE, stringsAsFactors = FALSE)
# check that plinkname is the same as popmap
if (!all(plinkname$V1 == popmap$SampleId)) {
    stop('Admixture file contains mismatching samples')
}

# load admixture dataframe ========
adprefix = paste0(indir, 'maf', mafcut, '/all50_pass_maf', mafcut)

dtlist = lapply(Klist, read_admixture, prefix = adprefix, nrep = 10, popmap = popmap, subpoporder = subpoporder)
forplotlist = lapply(dtlist, format_admaster)

# main --------
# formatting masterdt ========
forplot = dplyr::bind_rows(forplotlist)
# set factor levels for plotting (subpopulation ordered by latitude)
forplot$SubPopId = factor(forplot$SubPopId, levels = subpoporder)
forplot$SampleId = factor(forplot$SampleId, levels = get_sampleorder(popmap, subpoporder))

# plotting ========
pp <- ggplot(forplot, aes(x = SampleId, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set3") +
    facet_grid(K ~ run) +
    scale_y_continuous(labels = scales::percent,expand = c(0,0)) +
    labs(y = 'Ancestry Fraction') + 
    theme_pubclean() + 
    theme(axis.text.x = element_text(angle = 90, size = 6),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank()) 

# output files --------
ggsave(filename = paste0('Admixture_maf', mafcut, '_5runs_K26_', dataset, '_', ref, '_', today, '.pdf'), path = plotdir, 
       plot = pp, height = 8, width = 14)

saveRDS(forplot, file = paste0(outdir, 'Admixture_maf', mafcut, '_10runs_K26_', dataset, '_', ref, '_', today, '.rds'))
write.csv(forplot, file = paste0(outdir, 'Admixture_maf', mafcut, '_10runs_K26_', dataset, '_', ref, '_', today, '.csv'))
# cleanup --------
date()
closeAllConnections()
