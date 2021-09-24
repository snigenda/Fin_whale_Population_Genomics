# Title: Plot the site counts 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed May 20 15:08:49 2020

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

# no sorting!
get_topn <- function(df, topn) {
    # df keep 
    outdf = df[1:topn,]
    otherdf = df[(topn+1):nrow(df),] 
    otherdf = otherdf %>%  dplyr::summarise(across(where(is.numeric), sum))
    rownames(otherdf) = "Other" 
    outdf = base::rbind(outdf, otherdf) %>% 
        tibble::rownames_to_column(var = "SITE_FILTER")
    outdf$SITE_FILTER = factor(outdf$SITE_FILTER, levels = outdf$SITE_FILTER)
    return(outdf)
}

# def variables --------
# args = commandArgs(trailingOnly=TRUE)
# dataset = as.character(args[1])
# ref = as.character(args[2])
dataset = "all50"
ref = "Minke"
cutoff = 10 # the max percentage of reason to exclude 
today = format(Sys.Date(), "%Y%m%d")
datadate = "20201110"

# go to that dataset 
setwd(paste0("~/Lab/fin_whale/data/analyses/Summary_stats/", dataset))
plotdir = "./plots/"
outdir = "./summary_stats/"
dir.create(plotdir, recursive = T)
dir.create(outdir, recursive = T)

# load data --------
filtertally = data.frame()
for (ii in 1:get_contigs(ref)) {
    idx = stringr::str_pad(as.character(ii), width = 2, pad = "0")
    filterfile = paste0(ref, "/filter_stats_", datadate, "/", dataset, "_", ref, "_", idx,"_filter_stats_FILTER_tally.csv")
    fttemp = read.csv(filterfile, stringsAsFactors = F) 
    if (ii == 1) {
        filtertally = fttemp
    } else {
        filtertally = dplyr::full_join(filtertally, fttemp, by = "SITE_FILTER")
    }
}
# fill in zero and sort 
filtertally = filtertally %>% base::replace(is.na(.), 0)
filtertally$N_FILTER_all = rowSums(filtertally[,-1])
filtertally = filtertally %>% 
    dplyr::arrange(desc(N_FILTER_all))
# write the data out 
write.csv(filtertally, file = paste0(outdir, "filter_tally_", today, ".csv"))
sum(filtertally$N_FILTER_all)
# [1] 2324429748

# convert to percentage 
filterper = filtertally %>% 
    tibble::column_to_rownames(var = "SITE_FILTER") 
filterper = apply(filterper, 2, function(x) x / sum(x)) %>% as.data.frame()
write.csv(filterper, file = paste0(outdir, "filter_tally_per_", today, ".csv"))

# plot --------
filterpertop10 = get_topn(filterper, cutoff) %>% 
    reshape2::melt(id.vars = "SITE_FILTER")

pp1 <- ggplot(data = filterpertop10, aes(x = variable, y = value, fill = SITE_FILTER)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Contiglist", y = "Filter %", title = "Site filter distribution") +
    scale_x_discrete(labels = c(as.character(1:get_contigs(ref)), "all")) +
    scale_fill_brewer(palette = "Set3") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0("filter_tally_per_", today, ".pdf"), plot = pp1, path = plotdir, width = 12, height = 6)
