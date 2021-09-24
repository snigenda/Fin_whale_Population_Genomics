# Title: Plot the site counts 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Aug 16 22:02:11 2020

# preparation --------
rm(list = ls())
cat("\014")

library(dplyr)
library(reshape2)
library(ggplot2)

source("~/Lab/fin_whale/scripts/config/plotting_config.R")

# def functions --------
get_contigs <- function(ref) {
    if (ref == "Minke") {
        ncontig = 96
    } 
    if (ref == "Bryde") {
        ncontig = 23
    }
    return(ncontig)
}

# sort by the total order 
sort_gttally <- function(df) {
    # df keep 
    test = df %>% 
        tibble::column_to_rownames(var = "GT_FILTER")
    gtfilter = names(sort(rowSums(test), decreasing = T))
    df$GT_FILTER = factor(df$GT_FILTER, levels = gtfilter)
    return(df)
}

# def variables --------
# args = commandArgs(trailingOnly=TRUE)
# dataset = as.character(args[1])
# ref = as.character(args[2])
dataset = "all50"
ref = "Minke"
today = format(Sys.Date(), "%Y%m%d")
datadate = "20201110"

# go to that dataset 
setwd(paste0("~/Lab/fin_whale/data/analyses/Summary_stats/", dataset))
plotdir = "./plots/"
outdir = "./summary_stats/"
dir.create(plotdir, recursive = T)
dir.create(outdir, recursive = T)

# load data --------
gtfiltertally = data.frame()
for (ii in 1:get_contigs(ref)) {
    idx = stringr::str_pad(as.character(ii), width = 2, pad = "0")
    gtfilterfile = paste0(ref, "/filter_stats_", datadate, "/", dataset, "_", ref, "_", idx,"_filter_stats_GTFILTER_tally.csv")
    gtfttemp = read.csv(gtfilterfile, stringsAsFactors = F, 
                        colClasses = c("character", rep("numeric", times = 50))) %>% 
        dplyr::mutate(CONTIGID = idx) %>% 
        base::replace(is.na(.), 0)
    colnames(gtfttemp)[1] = "GT_FILTER"
    if (ii == 1) {
        gtfiltertally = gtfttemp
    } else {
        gtfiltertally = dplyr::bind_rows(gtfiltertally, gtfttemp)
    }
}

# check the rowSum should be the same as the total sites 
totalgt = colSums(apply(gtfiltertally[,2:51], 2, as.numeric))
unique(totalgt) # 2324429748 (the total genotype called available)

# write the raw data output 
write.csv(gtfiltertally, file = paste0(outdir, "gtfilter_tally_sepcontig_", today, ".csv"))

# group by filter types and discard the contigid info
gtfiltertallyall = gtfiltertally %>% 
    dplyr::select(-CONTIGID) %>% 
    dplyr::group_by(GT_FILTER) %>% 
    dplyr::summarise_all(sum)
write.csv(gtfiltertallyall, file = paste0(outdir, "gtfilter_tally_all_", today, ".csv"))

# plot --------
totalsites = unique(totalgt)
forplot = sort_gttally(gtfiltertallyall) %>% 
    reshape2::melt(id.vars = "GT_FILTER") %>% 
    dplyr::mutate(percent = value/totalsites)

pp1 <- ggplot(data = forplot, aes(x = variable, y = value, fill = GT_FILTER)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "SampleId", y = "Filter count", title = "Genotype filter distribution") +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("gtfilter_tally_all_", today, ".pdf"), plot = pp1, path = plotdir, width = 8, height = 6)

pp2 <- ggplot(data = forplot, aes(x = variable, y = percent, fill = GT_FILTER)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "SampleId", y = "Filter count", title = "Genotype filter distribution") +
    scale_fill_brewer(palette = "Set3") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0("gtfilter_tally_per_", today, ".pdf"), plot = pp2, path = plotdir, width = 8, height = 6)