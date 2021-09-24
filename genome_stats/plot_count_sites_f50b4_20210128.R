# Title: Plot the site counts for baleen whales as well
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Jan 28 14:56:29 2021

# Output: Tally of the total genotyped sites and basic statistics after variant filtration for f50b4 dataset
# Execution: Rscript --vanilla plot_count_sites_f50b4_20210128.R


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(reshape2)

source("~/Lab/finwhale_manuscript/scripts/config/plotting_config.R")

# def functions --------

# def variables --------
dataset = "f50b4"
nsample = 54
ref = "Minke"
today = format(Sys.Date(), "%Y%m%d")
datadate = "20210127" # the date the data generated

# go to that dataset
setwd('~/Lab/finwhale_manuscript/')
indir = './data/genome_stats/raw_data/'
outdir = './data/genome_stats/derived_data/'
dir.create(outdir)
sessionInfo()

# load data --------
summaryfiles = paste0('./data/genome_stats/raw_data/f50b4_count_sites_20210127/f50b4_Minke_', sprintf('%02d', seq(1,96)), '_sites_summary.txt')
sumdt_list = lapply(summaryfiles, read.delim, stringsAsFactors = F)
sumdt = dplyr::bind_rows(sumdt_list)

minke_contiglist = read.csv(file = "./scripts/config/minke_contig_summary.csv", row.names = 1, stringsAsFactors = F)
pop_map = read.csv(file = "./scripts/config/popmap_f50b4.csv", header = T, stringsAsFactors = F)

# main --------
sumdt = sumdt %>%
    dplyr::left_join(., y = pop_map, by = 'SampleId') %>%
    dplyr::mutate(HomAltRealCount = HomAltCount - HomAltAllCount,
                  ContigHet = HetCount/CalledCount)

# Basic stats --------
# more individuals, more sites failed variant filtration
PassCount = sum(unique(sumdt$PassCount))
# 890858824(all50_20201115) --> 880177286 (f50b4_20210127)
TotalCount = sum(unique(sumdt$TotalCount))
# 2324429748(all50_20201115) --> 2324429717 (f50b4_20210127)
TotalRefLen = sum(minke_contiglist$LN)
# 2324429847

# get the average heterozygosity
totaldt = sumdt %>%
    dplyr::group_by(SampleId, SpeciesId, PopId, SubPopId) %>%
    dplyr::summarise(TotalHomRef = sum(HomRefCount),
                     TotalHomAlt = sum(HomAltCount),
                     TotalHomAltAll = sum(HomAltAllCount),
                     TotalHet = sum(HetCount),
                     TotalCalled = sum(CalledCount),
                     TotalMissing = sum(MissingCount),
                     TotalPass = sum(PassCount),
                     TotalgVCF = sum(as.numeric(TotalCount)),
                     .groups = 'drop') %>%
    dplyr::mutate(GenomeHet = TotalHet/TotalCalled)
# These genomewide heterozygosity was used in Figure 2 plotting
# Also in the main text: "Compared with other marine mammals that have experienced different levels of population contractions, such as the diminutive vaquita in the Gulf of California35,36 (0.1 het/kb), abundant minke whale37 (0.6 het/kb) and endangered blue whale38 (2.1 het/kb)"

# output files --------
write.csv(totaldt, file = paste0(outdir, "f50b4_genomewide_heterozygosity_", today, ".csv"))
write.csv(sumdt, file = paste0(outdir, "f50b4_sites_summary_combined_", today, ".csv"))

# cleanup --------

