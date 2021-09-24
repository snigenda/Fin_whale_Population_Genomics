# Title: Subset the Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.bed to within the regions called
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 19:46:27 2020

# preparation --------
rm(list = ls())
cat("\014")

options(echo = TRUE)

# def functions --------

# def variables --------

# main -------
contiglist = read.csv(file = '/u/project/rwayne/snigenda/finwhale/scripts/config/minke_contig_summary.csv', row.names = 1, stringsAsFactors = FALSE)

bedfile = read.delim(file = '/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.bed', header = FALSE, stringsAsFactors = FALSE)

bedout = bedfile[bedfile$V1 %in% contiglist$SN, ]

# output --------
write.table(bedout, file = '/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/bedfiles/Minke.CDSstartstop.Coords.FromGTF.0Based.sorted_merged.inVCF.bed', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
date()
