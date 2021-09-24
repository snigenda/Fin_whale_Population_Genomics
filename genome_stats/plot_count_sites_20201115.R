# Title: Plot the site counts
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed May 20 15:08:49 2020

# Output: Tally of the total genotyped sites and basic statistics after variant filtration
# Execution: Rscript --vanilla plot_count_sites_20201115.R

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(reshape2)

source("~/Lab/finwhale_manuscript/scripts/config/plotting_config.R")

# def functions --------

# def variables --------
dataset = "all50"
nsample = 50
ref = "Minke"
ncontig = 96

today = format(Sys.Date(), "%Y%m%d")
datadate = "20201115" # the date the data generated

# go to that dataset
setwd('~/Lab/finwhale_manuscript/')
indir = './data/genome_stats/raw_data/'
outdir = './data/genome_stats/derived_data/'
dir.create(outdir)
sessionInfo()

# load data --------
summaryfiles = paste0('./data/genome_stats/raw_data/all50_count_sites_20201115/all50_Minke_', sprintf('%02d', seq(1,96)), '_sites_summary.txt')
sumdt_list = lapply(summaryfiles, read.delim, stringsAsFactors = F)
sumdt = dplyr::bind_rows(sumdt_list)

minke_contiglist = read.csv(file = "./scripts/config/minke_contig_summary.csv", row.names = 1, stringsAsFactors = F)

# main --------
sumdt = sumdt %>%
    dplyr::mutate(PopId = substr(SampleId, 1, 3),
                  SubPopId = ifelse(PopId == "ENP",
                                    substr(SampleId, 4, 5), PopId),
                  HomAltRealCount = HomAltCount - HomAltAllCount,
                  ContigHet = HetCount/CalledCount)

# Basic stats --------
PassCount = sum(unique(sumdt$PassCount))
# 890858824(all50_20201115)
TotalCount = sum(unique(sumdt$TotalCount))
# 2324429748(all50_20201115)
TotalRefLen = sum(minke_contiglist$LN)
# 2324429847
# Note the differences in TotalCount and TotalRefLen is mostly due to the errors during JointGenotyping such as too many alleles issues 

# get the average heterozygosity
totaldt = sumdt %>%
    dplyr::group_by(SampleId, PopId, SubPopId) %>%
    dplyr::summarise(TotalHomRef = sum(HomRefCount),
                     TotalHomAlt = sum(HomAltCount),
                     TotalHomAltAll = sum(HomAltAllCount),
                     TotalHet = sum(HetCount),
                     TotalCalled = sum(CalledCount),
                     TotalMissing = sum(MissingCount),
                     TotalPass = sum(PassCount),
                     TotalgVCF = sum(as.numeric(TotalCount)),
                     .groups = 'drop') %>%
    dplyr::mutate(GenomeHet = TotalHet/TotalCalled) # genomewide heterozygosity

# output files --------
write.csv(totaldt, file = paste0(outdir, "all50_genomewide_heterozygosity_", today, ".csv"))
write.csv(sumdt, file = paste0(outdir, "all50_sites_summary_combined_", today, ".csv"))

# cleanup --------

