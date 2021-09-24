# Title: List the vcfR files to use for generating matching results
# 1. Add a vcfR object that restrict to protein_coding and no warnings in snpEff
# 2. Remove invariant sites
# 3. Generate genotype table summary
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Apr 22 15:36:20 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(vcfR)
library(stringr)
library(dplyr)

source("~/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# def functions --------
# split ann
split_ann <- function(info) {
    dt = unlist(lapply(info, function(xx){
        yy = stringr::str_split(xx, pattern = ',')[[1]][1]
        return(yy)
    }))

    dt2 = reshape2::colsplit(dt, pattern = "\\|", names = annfields)
    return(dt2)
}

# extract pass
extract_pass <- function(ann) {
    print(table(ann$`ERRORS / WARNINGS / INFO`, useNA = 'always'))
    print(table(ann$Transcript_BioType, useNA = 'always'))
    passid = (ann$Transcript_BioType == 'protein_coding' & ann$`ERRORS / WARNINGS / INFO` == '')
    print(table(passid, useNA = 'always'))
    return(passid)
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
    print(dim(vcfgt))
    print(dim(outgt))
    return(outgt)
}

generate_summary <- function(gt, muttype) {
    # start with output df with sample name at each row and tally info at each column
    # CalledCount = "HomRefCount"+"HomAltCount"+"HetCount"
    # PassCount = "CalledCount"+"MissingCount"
    tallyinfo = c("SampleId","HomRefCount","HomAltCount","HetCount","CalledCount","MissingCount","PassCount")
    outdf = data.frame(matrix(nrow = ncol(gt), ncol = length(tallyinfo)))
    colnames(outdf) = tallyinfo
    outdf[,'SampleId'] = colnames(gt)
    # loop through each sample
    for(ii in 1:ncol(gt)) {
        temp = gt[,ii]
        # fix the na
        temp[is.na(temp)] = './.'
        # check if right sample
        if (outdf[ii,'SampleId']!=colnames(gt)[ii]) {
            stop('ERROR: SampleId mismatch!')
        }
        outdf[ii,'HomRefCount']=sum(temp == '0/0')
        outdf[ii, 'HomAltCount']=sum(temp == '1/1')
        outdf[ii, 'HetCount']=sum(temp == '0/1')
        outdf[ii, 'CalledCount']=sum(temp != './.')
        outdf[ii, 'MissingCount']=sum(temp == './.')
        # sanity check
        if (outdf[ii, 'CalledCount'] != outdf[ii,'HomRefCount'] + outdf[ii, 'HomAltCount'] + outdf[ii, 'HetCount']) {
            stop('ERROR: CalledCount mismatch!')
        }
        # sanity check
        if (nrow(gt) != outdf[ii,'CalledCount'] + outdf[ii, 'MissingCount']) {
            stop('ERROR: PassCount mismatch!')
        }
    }
    outdf[,'PassCount'] = nrow(gt)
    outdf = outdf %>%
        dplyr::mutate(HomRefPer = HomRefCount/CalledCount,
                      HomAltPer = HomAltCount/CalledCount,
                      HetPer = HetCount/CalledCount,
                      RefAlleleCount = 2*HomRefCount + HetCount,
                      AltAlleleCount = 2*HomAltCount + HetCount,
                      CalledAlleleCount = 2*CalledCount,
                      RefAllelePer = RefAlleleCount/CalledAlleleCount,
                      AltAllelePer = AltAlleleCount/CalledAlleleCount,
                      PopId = substr(SampleId, 1, 3),
                      SubPopId = ifelse(PopId == "ENP",
                                        substr(SampleId, 4, 5), PopId),
                      MutType = muttype)
    # sanity check: CalledAlleleCount = RefAlleleCount + AltAlleleCount
    if (!(all(outdf$CalledAlleleCount == outdf$RefAlleleCount + outdf$AltAlleleCount))) {
        stop('ERROR: Allele Count mismatch!')
    }
    return(outdf)
}


# wrapper of generate_summary
wrapper_summary <- function(gtlist, gtprefix = 'SUMtable_Seg_PASSm6') {
    muttypes = names(gtlist)
    tallylist = vector(mode = 'list', length = length(gtlist))
    for (ii in 1:length(gtlist)) {
        mut = muttypes[ii]
        cdstype = cdstypelist[ii]
        gt = gtlist[[ii]]
        mytally = generate_summary(gt = gt, muttype = mut)
        # output tally file
        write.csv(mytally, file = paste0(outdir3,
                                         paste(mut, dataset, ref, cdstype, gtprefix, today, sep = '_'),
                                         '.csv'))
        tallylist[[ii]] = mytally
    }
    sumtally = dplyr::bind_rows(tallylist)
    return(sumtally)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50_snpEff_matching'
ref = 'Minke'
cdstypelist = c('filteredvcf','ALLregions','ALLregions','filteredvcf')
prefixlist = c('HIGH','MODERATE','LOW','MODIFIER')
projectpth = '/Users/linmeixi/google_drive/finwhale/analyses/'
vcfRpathes = c('snpEff_impact/all50/Minke/derive_data/vcfR/HIGH_all50_Minke_filteredvcf_PASSm6_20210416.rds',
               'DelVar_vcfR/all50_snpEff_impactCDS/Minke/derive_data/vcfR/MODERATE_all50_snpEff_impactCDS_Minke_ALLregions_PASSm6_20210415.rds',
               'DelVar_vcfR/all50_snpEff_impactCDS/Minke/derive_data/vcfR/LOW_all50_snpEff_impactCDS_Minke_ALLregions_PASSm6_20210415.rds')
annfields = c('Allele','Annotation','Annotation_Impact','Gene_Name','Gene_ID','Feature_Type','Feature_ID','Transcript_BioType','Rank','HGVS.c','HGVS.p','cDNA.pos / cDNA.length','CDS.pos / CDS.length','AA.pos / AA.length','Distance','ERRORS / WARNINGS / INFO')

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR', dataset, ref, sep = '/')
dir.create(workdir, recursive = TRUE)
setwd(workdir)

outdir1 = './derive_data/vcfR/'
outdir2 = './derive_data/gttable/'
outdir3 = './derive_data/sum_table/'
logdir = "./logs/"

dir.create(outdir1, recursive = TRUE)
dir.create(outdir2, recursive = TRUE)
dir.create(outdir3, recursive = TRUE)
dir.create(logdir, recursive = FALSE)

today = format(Sys.Date(), "%Y%m%d")
# sink(file = paste0(logdir, 'step1_listvcfR_', paste(dataset, ref, today, sep = '_'), ".log"))

date()
sessionInfo()

# load data --------
vcflist = lapply(vcfRpathes, function(xx) {
    yy <- readRDS(file = paste0(projectpth, xx))
    print(dim(yy))
    return(yy)
})

names(vcflist) = prefixlist[1:length(vcflist)]

# main --------
annlist <- lapply(vcflist, function(vcf) {
    # get INFO
    info = vcfR::extract.info(vcf, element = 'ANN')
    ann = split_ann(info)
    return(ann)
})

# get sites that are protein coding and has no warning ========
passidlist <- lapply(annlist, function(ann) {
    # get passid
    passid = extract_pass(ann)
    return(passid)
})

# subsetvcfs ========
vcfCLlist = vector(mode = 'list', length = length(vcflist))
for (ii in 1:length(vcflist)) {
    vcf = vcflist[[ii]]
    passid = passidlist[[ii]]
    vcf2 = vcf[passid,]
    print(dim(vcf2))
    vcfCLlist[[ii]] = vcf2
}

names(vcfCLlist) = names(vcflist)

# get genotype tables ========
gtlist <- lapply(vcflist, function(vcf) {
    gt = vcfR::extract.gt(vcf)
    gt_seg = remove_invar_sites(gt)
    return(gt_seg)
})

gtCLlist <- lapply(vcfCLlist, function(vcf) {
    gt = vcfR::extract.gt(vcf)
    gt_seg = remove_invar_sites(gt)
    return(gt_seg)
})

# get genotype tally ========
tallydt = wrapper_summary(gtlist, gtprefix = 'SUMtable_Seg_PASSm6')
tallyCLdt = wrapper_summary(gtCLlist, gtprefix = 'SUMtable_Seg_PASSm6CL')

# output files --------
# save file for the PASSm6 dataset ========
# save the vcflist  with passing sites and removed individuals
saveRDS(object = vcflist, file = paste0(outdir1, paste('HML', dataset, ref,'PASSm6', today, sep = '_'), '.rds'))
# 'PASSm6' and remove all invariant sites accounting for missing sites
saveRDS(object = gtlist, file = paste0(outdir2,  paste('HML', dataset, ref,'GTtable_Seg_PASSm6', today, sep = '_'), '.rds'))
# the summary dataframe
saveRDS(object = tallydt, file = paste0(outdir3,  paste('HML', dataset, ref,'SUMtable_Seg_PASSm6', today, sep = '_'), '.rds'))

# save file for the CLeaned PASSm6 dataset ========
# save the vcflist  with passing sites and removed individuals and sites in protein coding regions and sites without warnings in snpEff annotations
saveRDS(object = vcfCLlist, file = paste0(outdir1, paste('HML', dataset, ref,'PASSm6CL', today, sep = '_'), '.rds'))
# 'PASSm6' and remove all invariant sites accounting for missing sites
saveRDS(object = gtCLlist, file = paste0(outdir2,  paste('HML', dataset, ref,'GTtable_Seg_PASSm6CL', today, sep = '_'), '.rds'))
# the summary dataframe
saveRDS(object = tallyCLdt, file = paste0(outdir3,  paste('HML', dataset, ref,'SUMtable_Seg_PASSm6CL', today, sep = '_'), '.rds'))

# cleanup --------
date()
# sink()
closeAllConnections()

