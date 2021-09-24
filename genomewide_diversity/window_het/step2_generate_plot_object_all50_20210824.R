# Title: Window-based Heterozygosity plot objects
# Author: Jacqueline Robinson, Paulina Nunez Valencia (pnunez@lcg.unam.mx), Sergio Nigenda and Meixi Lin (meixilin@ucla.edu)
# Date: Wed Mar 24 09:53:52 2021

# Output: formatted `rds` object which is flexible for plotting
# Execution: Rscript --vanilla step2_generate_plot_object_all50_20210325.R

# preparation --------
rm(list = ls())
cat("\014")

options(echo = TRUE)

library(plyr)
library(dplyr)

source("~/Lab/finwhale_manuscript/scripts/config/plotting_config.R")

# def functions --------
# get seq along coord for plotwinhet_merge (NOT MATCHING THE REAL COORDINATES BUT CONSIDERS THE END OF EACH WINDOW IN SCAFFOLDS)
get_merge_start <- function(dt) {
    merge_start = dt[1,'window_start'] # initiate
    for (ii in 2:nrow(dt)) {
        next_start = merge_start[length(merge_start)] + dt[ii-1, 'winlen']
        merge_start=c(merge_start, next_start)
    }
    dt = cbind(dt, merge_start)
    return(dt)
}

# format the temphet file
format_temphet <- function(temphet,sample,subpop) {
    colnames(temphet) = c("chrom","window_start","window_end","sites_total","calls", "hets")
    forplot <- temphet %>%
        dplyr::mutate(hetpkb = 1000*hets/calls,
                      winlen = window_end - window_start + 1)
    forplot <- get_merge_start(forplot) %>%
        dplyr::mutate(sample = sample,
                      subpop = subpop)
    return(forplot)
}

# def variables --------
dataset = "all50"
nsample = 50
ref = "Minke"
today = format(Sys.Date(), "%Y%m%d")
datadate = "20201206" # the date the data generated

# window het settings
winsize = 1e+6
stepsize = 1e+6
misscutoff = 0.2 # minimum called/winsize (same as isle royale paper)

# go to that dataset
setwd('~/Lab/finwhale_manuscript/')
indir = './data/window_het/raw_data/all50_window_het_20201206/'
outdir = './data/window_het/derived_data/'
dir.create(outdir)
sessionInfo()

# load data --------
contiglist = read.csv(file = "./scripts/config/minke_contig_summary.csv", row.names = 1, stringsAsFactors = F)
pop_map = read.csv(file = "./scripts/config/popmap_all50.csv", header = T, stringsAsFactors = F)

# load the heterozygosity summary
hetfiles = paste0('./data/window_het/raw_data/all50_window_het_20201206/JointCalls_all50_08_B_VariantFiltration_', sprintf('%02d', seq(01,96)), '_het_1000000win_1000000step.txt')
allhet <- plyr::ldply(hetfiles, read.table, stringsAsFactors = FALSE, header=TRUE, sep="\t")

# the regions with too much missing sites are the same across individuals
# remove the regions with too much missing or the region at the end of a scaffold that are way too short compared with 1mb window (should be 2284)
cleanhet <- allhet[which(allhet[,'sites_total']>=(misscutoff*winsize)),] 

# main -------
# first get data for plotting ========
plotdt = vector(mode = "list")
for(ii in 1:nsample) {
    sample <- pop_map[ii,]
    barcol <- loccolors[sample$SubPopId]
    temphet <- cleanhet %>% dplyr::select(chrom,window_start,window_end,sites_total,ends_with(sample$SampleId))
    forplot <- format_temphet(temphet, sample = sample$SampleId, subpop = sample$SubPopId)
    plotdt[[ii]] = forplot
}

plotdt = dplyr::bind_rows(plotdt)

# output data --------
write.csv(allhet, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_all50_', today, '.csv'))
saveRDS(plotdt, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_', today, '.rds'))
write.csv(plotdt, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_', today, '.csv'))

# cleanup --------
closeAllConnections()
date()
