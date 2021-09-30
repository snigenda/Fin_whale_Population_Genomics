# Title: Use the ape neighbor joining tree to construct phylogeny (For use in hoffman2 cluster for bootstrapping, no ggplotting part)
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Mon Mar  1 02:15:28 2021
# Example: Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/f50b4/step3_ApePhylogeny_f50b4_20210301.R '/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf05.gds' 'f50b4_pass_bialleic_all_LDPruned_maf05'

# IMPORTANT NOTE: The input rows (sample names) for ape::nj need to be named with values other than 1,2,3... The dist.gene function will confuse the dimmension of genotype
# IMPORTANT NOTE: Here the input matrix should be on the dosage of alternative alleles (when you use the dosage of reference allele, tree looks differently.)
# IMPORTANT NOTE: The ggtree is pretty outdated in bioconda. https://github.com/YuLab-SMU/ggtree/issues/91

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(gdsfmt)
library(SNPRelate)
# library(SeqArray)
library(ape)
library(dplyr)

# source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R')
source('/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/plotting_config.R')

# def functions --------
# convert gtdt to the input for ape::dist.gene
get_njdt <- function(genofile) {
    gtdt <- read.gdsn(index.gdsn(genofile, "genotype"))
    sampleid <- read.gdsn(index.gdsn(genofile, "sample.id"))
    dimnames(gtdt)[1] <- list(sampleid)
    gtdt = 2-gtdt # convert the dosage to alternative alleles
    gtdt[gtdt==-1] <- NA
    return(gtdt)
}

# get maf
get_maf <- function(gdsfile) {
    maf = stringr::str_split(gdsfile, pattern = '_')[[1]]
    maf = maf[length(maf)]
    maf = stringr::str_split(maf, pattern = '.gds')[[1]][1]
    return(maf)
}

# get neighbor-joining trees along with the bootstrap trees
get_trees <- function(njdt, trmethod, nbt = 1000) {
    ### Create trees
    if(trmethod == "percentage") {
        #  The rows of the data matrix represent the individuals, and the columns the loci.
        distmat <- ape::dist.gene(njdt, method = trmethod, pairwise.deletion = TRUE, variance = FALSE)
        f <- function(x) {ape::nj(ape::dist.gene(x, method = trmethod, pairwise.deletion = TRUE, variance = FALSE))}
    }
    else {
        distmat <- ape::dist.gene(njdt, method = trmethod, pairwise.deletion = FALSE, variance = FALSE)
        f <- function(x) {ape::nj(ape::dist.gene(x, method = trmethod, pairwise.deletion = FALSE, variance = FALSE))}
    }

    tree_nj <- ape::nj(distmat)
    btrees <- ape::boot.phylo(tree_nj, njdt, f, B = nbt, quiet = F, trees = T)
    tree_nj$node.label <- btrees$BP
    output = list(tree_nj, btrees$trees)

    return(output)
}

# def variables --------
args <- commandArgs(trailingOnly=TRUE)
workdir <- as.character(args[1]) # the working directory (should be the same as before and after)
gdsfile <- as.character(args[2]) # the name of gdsfile to analyze
outprefix <- as.character(args[3]) # the output prefix

today = format(Sys.Date(), "%Y%m%d")
outgrp = 'EubGla01'
nbt = 1000
maf = get_maf(gdsfile)
outdir = './derive_phylogeny/'

setwd(workdir)
sessionInfo()
dir.create(outdir)

# load data --------
snpgdsSummary(gdsfile) # print summary
genofile <- SNPRelate::snpgdsOpen(gdsfile, readonly = TRUE)
# extract the genotypes
njdt <- get_njdt(genofile)

# get the population names
# popmap = read.csv(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/baleen_popid.csv', stringsAsFactor = FALSE)
popmap = read.csv(file = '/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/baleen_popid.csv', stringsAsFactor = FALSE)

# main --------
# get the ape trees ========
# first object is the tree, second object is the boostrapped trees
treespair <- get_trees(njdt, trmethod = 'pairwise', nbt = nbt)
treesperc <- get_trees(njdt, trmethod = 'percentage', nbt = nbt)

# output files --------
saveRDS(treespair, file = paste0(outdir, outprefix, '_treespair_', today, '.rds'))
saveRDS(treesperc, file = paste0(outdir, outprefix, '_treesperc_', today, '.rds'))

# cleanup --------
closefn.gds(genofile)
date()
closeAllConnections()
