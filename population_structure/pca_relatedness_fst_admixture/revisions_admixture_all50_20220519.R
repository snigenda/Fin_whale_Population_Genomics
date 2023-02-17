# Title: Plot admixture output for all scenarios in all50 populations
# Plan to be included in supp figures
# Author: Meixi Lin (meixilin@ucla.edu)
# Thu May 19 10:47:01 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

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
        if (!file.exists(qfile[ii])){
            print(qfile[ii])
            stop('File not exist')
        }
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
subpoporder = c("AK", "BC", "WA", "OR", "CA") # order for subpopulations

gdsfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds'

# set derived data outside of the 'Minke' folder to make hoffman2 folder unique
workdir = '/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/PopStructure/all50/'
setwd(workdir)
indir = './Minke/Admixture_20210318/maf10/'
outdir = './derive_data/revisions_Admixture_all50_maf10/'
plotdir = './plots/revisions_Admixture_all50_maf10/'

dir.create(outdir)
dir.create(plotdir)

# load data --------
popmap = read.csv(file = "/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/popmap_all50.csv", stringsAsFactors = F)

plinkname = read.table(file = './Minke/JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10_SA_mrF.nosex', header = FALSE, stringsAsFactors = FALSE)

# check that plinkname is the same as popmap
if (!all(plinkname$V1 == popmap$SampleId)) {
    stop('Admixture file contains mismatching samples')
}

# load admixture dataframe ========
adprefix = paste0(indir, 'all50_pass_maf10')

dtlist = lapply(Klist, read_admixture, prefix = adprefix, nrep = 10, popmap = popmap, subpoporder = subpoporder)
forplotlist = lapply(dtlist, format_admaster, keeprun = 1:10)

# load CV data ========
# read cross validation files
adcv = read.csv(file = paste0(indir, 'Admixture_CV_LLsummary_maf10.csv'), stringsAsFactors = FALSE)

# main --------
# formatting masterdt ========
forplotQ = dplyr::bind_rows(forplotlist)
# set factor levels for plotting (subpopulation ordered by latitude)
forplotQ$SubPopId = factor(forplotQ$SubPopId, levels = subpoporder)
forplotQ$SampleId = factor(forplotQ$SampleId, levels = get_sampleorder(popmap, subpoporder))

# plotting ========
pp <- ggplot(forplotQ, aes(x = SampleId, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_vline(xintercept = 12.5, linetype = 'dotted') + 
    geom_vline(xintercept = 21.5, linetype = 'dotted') +
    geom_vline(xintercept = 30.5, linetype = 'dotted') +
    scale_fill_brewer(palette = "Set3") +
    facet_grid(run ~ K) +
    scale_y_continuous(labels = scales::percent,expand = c(0,0)) +
    labs(y = 'Ancestry Fraction') +
    theme_pubclean() +
    theme(axis.ticks.x = element_line(size = 0.3),
          axis.text.x = element_text(angle = 90, size = 4),
          axis.text.y = element_text(size = 7),
          strip.text.y = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank())

# Not combining plots since it was a pain ========
ggsave(filename = 'revisions_AdmixtureQ_all50.pdf', path = plotdir, plot = pp, height = 8, width = 12)

# output files --------
write.csv(forplotQ, file = paste0(outdir, 'revisions_admixture_all50_20220519.csv'))
save.image(file = paste0(outdir, 'revisions_admixture_all50_20220519.RData'))

# cleanup --------
date()
closeAllConnections()
