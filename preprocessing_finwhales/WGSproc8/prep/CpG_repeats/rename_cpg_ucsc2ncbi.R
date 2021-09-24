# Title: convert CpG islands from ucsc to ncbi format 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Mar 15 20:22:39 2020

# Usage: rename_cpg_ucsc2ncbi.R working_dir cpgIslandExtUnmasked.txt ucscToRefSeq.txt ref.dict out.name

# preparation --------
rm(list = ls())
cat("\014")
library(dplyr)

# def functions --------
	
# def variables --------
args = commandArgs(trailingOnly=TRUE)
mywd = as.character(args[1])
cpgisland = as.character(args[2])
ucsc2ncbi = as.character(args[3])
ref.dict = as.character(args[4])
outname = as.character(args[5])

cpgfields = c("bin", "chrom", "chromStart.cpg","chromEnd.cpg", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
ncbi2fields = c("chrom", "chromStart","chromEnd", "NCBIchrom")

# note the variables 
# ref.dict = "/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.dict"
# outname = "CpG_ExtUnmasked.bed"

# main -------
setwd(mywd)

# load the cpgisland 
cpg = read.delim(file = cpgisland, stringsAsFactors = F, header = F)
colnames(cpg) = cpgfields

# load the namechanges 
ncbi2 = read.delim(file = ucsc2ncbi, stringsAsFactors = F, header = F)
colnames(ncbi2) = ncbi2fields

# output joined bedfiles 
output = dplyr::left_join(cpg, ncbi2, by = "chrom") %>%
	dplyr::select(NCBIchrom, chromStart.cpg, chromEnd.cpg) %>%
	dplyr::tbl_df() %>% 
	dplyr::arrange(NCBIchrom)

# check if it encompassed all the contigs 
dict = read.delim(file = ref.dict, header = F, skip = 1, stringsAsFactors = F) %>%
    dplyr::tbl_df() %>% 
    dplyr::select(V2, V3, V4) %>%
    dplyr::mutate(SN = stringr::str_sub(V2, start = 4),
                  LN = as.integer(stringr::str_sub(V3, start = 4)))


# output --------
if (any(output$NCBIchrom %in% dict$SN) != TRUE) {
	print ("Error!! not all contigs were included.")
} else {
	write.table(output, file = outname, sep = "\t", quote = F, col.names = F, row.names = F)
}