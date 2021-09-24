# Title: Get genomewide coordinates for minke genome 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Sep 14 14:36:02 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/fin_whale/scripts/config/")

library(dplyr)
# def functions --------

# def variables --------

# load data --------
contiglist = read.csv("minke_contig_summary.csv", row.names = 1, stringsAsFactors = F, 
                      colClasses = c("character","character","character","character","character","integer", "character", "integer", "integer"))

# main --------
genomewide_coord = 1 # initiate
for (ii in 2:nrow(contiglist)) {
    temp = genomewide_coord[length(genomewide_coord)] + contiglist[ii-1, 'LN']
    genomewide_coord=c(genomewide_coord, temp)
}

# check for total length 
(genomewide_coord[nrow(contiglist)] + contiglist[nrow(contiglist), 'LN'] - 1) == sum(contiglist$LN)

# write the new contiglist 
contiglist = cbind(contiglist, genomewide_coord)
plot(x = contiglist$binid, y = contiglist$genomewide_coord, type = "b")
write.csv(contiglist, file = "minke_contig_summary.csv")

# cleanup --------
