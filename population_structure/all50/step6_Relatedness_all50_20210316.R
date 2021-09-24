# Title: Analyze Relatedness using the new LDPruned, output tables and figures
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Fri Mar 19 00:42:32 2021
# execute with: 
# source(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/all50/step6_Relatedness_all50_20210316.R', echo = TRUE, max.deparse.length = 1000)

# preparation --------
rm(list = ls())
cat("\014")

library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggpubr)
library(dplyr)

source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R')

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

setwd('/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/all50')
outdir = './derive_data/'
plotdir = './plots/'
logdir = './logs/'

dir.create(outdir)
dir.create(plotdir)
dir.create(logdir)

sink(file = paste0(logdir, 'Relatedness_all50_Minke_', today, '.log'))
dataset = 'all50'
ref = 'Minke'
mafcut = '10' # mafcutoff across the population used in the input gdsfile (NOT the maf setting used in the other samples)
missing = 0.2
gdsfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds'

sessionInfo()
getwd()

# def function --------
plot_kinship <- function(ibd.coeff) {
    pp <- ggplot(ibd.coeff, aes(x=ID1, y=ID2, fill=kinship)) + 
        geom_tile() + 
        theme_light() +
        scale_fill_gradient(low = '#1F618D', high = '#2ECC71') + 
        theme(axis.text.x=element_text(angle=90))
    return(pp)
}

# load data --------
genofile = SNPRelate::snpgdsOpen(filename = paste0(ref, '/', gdsfile))

#List of sample ids
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# sample designations
popmap = read.csv(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/popmap_all50.csv', stringsAsFactors = FALSE) 

# main --------
# estimate IBD in ENP ========
ENPids = popmap[popmap$PopId == 'ENP','SampleId']
# Estimating IBD Using PLINK method of moments (MoM)
# Missing rate and MAF set as the same in the tutorial (http://corearray.sourceforge.net/tutorials/SNPRelate/#f_st-estimation) and same as Paulina's previous settings
ENPibd <- snpgdsIBDMoM(genofile, sample.id = ENPids, verbose = TRUE,
                       maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE)
# k0	the probability of sharing ZERO alleles
# k1	the probability of sharing ONE alleles
# kinship	kinship coefficient
ENPibd.coeff <- snpgdsIBDSelection(ENPibd)
length(unique(ENPibd.coeff$ID1))
length(unique(ENPibd.coeff$ID2))

# get mean kinship
mean(ENPibd.coeff$kinship)

# estimate IBD in GOC ========
GOCids = popmap[popmap$PopId == 'GOC','SampleId']
# Estimating IBD Using PLINK method of moments (MoM)
# Missing rate and MAF set as the same in the tutorial (http://corearray.sourceforge.net/tutorials/SNPRelate/#f_st-estimation)
# Changing missing rate settings did not change the pattern overall 
GOCibd <- snpgdsIBDMoM(genofile, sample.id = GOCids, verbose = TRUE,
                       maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE)
# k0	the probability of sharing ZERO alleles
# k1	the probability of sharing ONE alleles
# kinship	kinship coefficient
GOCibd.coeff <- snpgdsIBDSelection(GOCibd)
length(unique(GOCibd.coeff$ID1))
length(unique(GOCibd.coeff$ID2))

# get mean kinship
mean(GOCibd.coeff$kinship)

# plot kinship matrix ========
ppENP <- plot_kinship(ibd.coeff = ENPibd.coeff)
ppGOC <- plot_kinship(ibd.coeff = GOCibd.coeff)

pp <- ggarrange(ppENP, ppGOC, nrow = 2, labels = c('A', 'B'))

ggsave(filename = paste0('Relatedness_maf', mafcut, '_SUBmaf05miss05_', dataset, '_', ref, '_', today, '.pdf'), path = plotdir, 
       plot = pp, height = 10, width = 6)

# output files --------
saveRDS(object = ENPibd, file = paste0(outdir, 'ENPibd_maf', mafcut, '_SUBmaf05miss05_', dataset, '_', ref, '_', today, '.rds'))
saveRDS(object = GOCibd, file = paste0(outdir, 'GOCibd_maf', mafcut, '_SUBmaf05miss05_', dataset, '_', ref, '_', today, '.rds'))
write.csv(ENPibd.coeff, file = paste0(outdir, 'ENPibdcoeff_maf', mafcut, '_SUBmaf05miss05_', dataset, '_', ref, '_', today, '.csv'))
write.csv(GOCibd.coeff, file = paste0(outdir, 'GOCibdcoeff_maf', mafcut, '_SUBmaf05miss05_', dataset, '_', ref, '_', today, '.csv'))

# cleanup --------
closefn.gds(genofile)
date()
sink()
closeAllConnections()

# load the files for additional validations --------
options(echo = TRUE)
ENPibd.coeff = read.csv(file = 'derive_data/ENPibdcoeff_maf10_SUBmaf05miss05_all50_Minke_20210323.csv', stringsAsFactors = FALSE)
summary(ENPibd.coeff$kinship)
ENPkin = ENPibd.coeff$kinship[ENPibd.coeff$kinship > 0]
summary(ENPkin)

GOCibd.coeff = read.csv(file = 'derive_data/GOCibdcoeff_maf10_SUBmaf05miss05_all50_Minke_20210323.csv', stringsAsFactors = FALSE)
summary(GOCibd.coeff$kinship)
GOCkin = GOCibd.coeff$kinship[GOCibd.coeff$kinship > 0]
summary(GOCkin)

boxplot(ENPkin, GOCkin)
boxplot(ENPibd.coeff$kinship, GOCibd.coeff$kinship)
# not removing zero
wilcox.test(ENPibd.coeff$kinship, GOCibd.coeff$kinship)
# Wilcoxon rank sum test with continuity correction
# 
# data:  ENPibd.coeff$kinship and GOCibd.coeff$kinship
# W = 1570, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

# removing zero values
wilcox.test(ENPkin, GOCkin)
# Wilcoxon rank sum test with continuity correction
# 
# data:  ENPkin and GOCkin
# W = 1322, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

