# Title: Read in the vcf files and store vcfR objects 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Feb 17 15:00:01 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(vcfR)
library(stringr)
library(dplyr)

source("~/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# def functions --------
# get good sample indices by providing bad sample names
# WARNING: Note that the first column of the ‘gt’ slot in `vcfR-vcf` contains the format information for all of the subsequent columns. This means you will typically want to include the first column when subsetting `vcfR-vcf`.
goodsample_index <- function(vcf, badnames) {
    samplevcf = colnames(vcfR::extract.gt(vcf))
    badids = unlist(lapply(badnames, function(xx) {which(samplevcf == xx)}))
    goodids = setdiff(1:length(samplevcf), badids)+1
    # sanity check
    if (length(goodids) != (length(samplevcf) - length(badids))) {
        stop('Wrong id lists!')
    }
    return(goodids)
}

# filter based on the FILTER=PASS
# can also discard individuals by `badsamples`, input NA if you don't want to subset
subset_vcf_filter_sample <- function(vcf,myfilter='PASS',badsamples=character(0)) {
    vcffilter=getFILTER(vcf)
    passid <- (vcffilter==myfilter)
    if (length(badsamples)==0) {
        outvcf <- vcf[passid,]
    } else {
        goodsample <- goodsample_index(vcf, badsamples)
        # include the first format field as well
        outvcf <- vcf[passid,c(1,goodsample)]
    }
    print(outvcf)
    return(outvcf)
}

# get invariant sites
count_invar_sites <- function(vcfgt) {
    invarcount11 = 0
    invarcount00 = 0
    for (ii in 1:nrow(vcfgt)) {
        temp=vcfgt[ii,]
        temp=temp[!is.na(temp)]
        if(all(temp == '1/1')) {
            # print(temp)
            invarcount11 = invarcount11 + 1
        } else {
            if (all(temp == '0/0')) {
                invarcount00 = invarcount00 + 1
            }
        }
    }
    invarcount = invarcount00 + invarcount11
    output = list(invarcount, invarcount00, invarcount11)
    return(output)
}

count_invar_sites2 <- function(vcfgt) {
    invarcount11 = 0
    invarcount00 = 0
    for (ii in 1:nrow(vcfgt)) {
        temp=vcfgt[ii,]
        temp[is.na(temp)] = './.'
        if(all(temp == '1/1')) {
            # print(temp)
            invarcount11 = invarcount11 + 1
        } else {
            if (all(temp == '0/0')) {
                invarcount00 = invarcount00 + 1
            }
        }
    }
    invarcount = invarcount00 + invarcount11
    output = list(invarcount, invarcount00, invarcount11)
    return(output)
}

# get invariant sites tally
get_invar_tally <- function(gt) {
    tally1=count_invar_sites(gt)
    tally2=count_invar_sites2(gt)
    print(paste('Accounting for ./.: Total =', tally1[[1]], 'Total 0/0 =', tally1[[2]], 'Total 1/1 =', tally1[[3]]))
    print(paste('NOT Accounting for ./.: Total =', tally2[[1]], 'Total 0/0 =', tally2[[2]], 'Total 1/1 =', tally2[[3]]))
}

# remove invariant sites (accounting for ./.)
remove_invar_sites <- function(vcfgt) {
    invarid = integer(0)
    for (ii in 1:nrow(vcfgt)) {
        temp=vcfgt[ii,]
        temp=temp[!is.na(temp)]
        if(all(temp == '1/1') | all(temp == '0/0')) {
            invarid = c(invarid, ii)
        }
    }
    passid = setdiff(1:nrow(vcfgt), invarid)
    outgt = vcfgt[passid,]
    return(outgt)
}

# def variables --------
args = commandArgs(trailingOnly=TRUE)
vcfprefix = as.character(args[1])
vcffile = as.character(args[2])
# vcfprefix = 'syn'
# vcffile = 'JointCalls_all50_filterpassmiss_syn_ALLregions_all_SIFT.vcf.gz'

dataset = 'all50'
ref = 'Minke'
cdstype = 'ALLregions'
indir = paste('/Users/linmeixi/google_drive/finwhale/analyses/get_ALLregions_CDS/', dataset, ref, sep = '/')
workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR/', dataset, ref, sep = '/')
dir.create(workdir, recursive = T)

setwd(workdir)
outdir1 = "./derive_data/vcfR/"
outdir2 = "./derive_data/gttable/"
logdir = "./logs/"
dir.create(outdir1, recursive = T)
dir.create(outdir2, recursive = T)
dir.create(logdir, recursive = F)

today = format(Sys.Date(), "%Y%m%d")
# sink(file = paste0(logdir, 'step1_createvcfR_', paste(vcfprefix, dataset, ref, cdstype, today, sep = '_'), ".log"))
date()
sessionInfo()

# load data --------
rawvcf=vcfR::read.vcfR(file = paste0(indir, "/", vcffile))

# main --------
# make a subset on only passing sites AND discard 6 individuals ========
# only pass
rawvcfP=subset_vcf_filter_sample(vcf = rawvcf)
# only pass + 44 passing samples
print(discardind)
rawvcfPS = subset_vcf_filter_sample(vcf = rawvcf, badsamples = discardind)

# get the genotype tables ========
gtP = vcfR::extract.gt(rawvcfP)
gtPS = vcfR::extract.gt(rawvcfPS)

# get the invariant sites tally ========
print('INFO: Getting invariant sites tally for only PASS ...')
get_invar_tally(gtP)
print('INFO: Getting invariant sites tally for only PASS + remove 6 individuals ...')
get_invar_tally(gtPS)

# remove all the invariant sites for the PASSm6 dataset ========
gtPS.seg = remove_invar_sites(gtPS)

print(dim(gtPS))
print(dim(gtPS.seg))

# save file ========
# for now only outputing the 'PASSm6' option
saveRDS(object = rawvcfPS, file = paste0(outdir1, paste(vcfprefix, dataset, ref, cdstype, 'PASSm6', today, sep = '_'), '.rds'))
# 'PASSm6' and remove all invariant sites accounting for missing sites
saveRDS(object = gtPS.seg, file = paste0(outdir2, paste(vcfprefix, dataset, ref, cdstype, 'GTtable_Seg_PASSm6', today, sep = '_'), '.rds'))

# cleanup --------
date()
# sink()
closeAllConnections()
