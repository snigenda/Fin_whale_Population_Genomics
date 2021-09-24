# Title: Prepare read statistics for variant filtering using files from
# prep_filter_vcf_gtdp_fast.sh
# Using all contigs information
# Author: Meixi Lin (meixilin@ucla.edu)
# Date:

# Usage:
# Rscript --vanilla parse_vcf_gtdp.R <FILEDIR>
# These files should exist: Minke_chr01_countgt.tsv; Minke_chr01_sumDP.tsv

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
library(plyr)
library(dplyr)

# def functions --------
# read files
namednum2df <- function(namednum, mycolnames) {
	data = as.data.frame(namednum) %>%
		tibble::rownames_to_column()
	colnames(data) = mycolnames
	return(data)
}

# def variables --------
# get input
args = commandArgs(trailingOnly=TRUE)
mywd = as.character(args[1])

setwd(mywd)

# preset variables
dataset = "f50b4"
ref = "Minke"
datadate = "20210124"
lower = 1/3
upper = 2.5 # upper limit of files
nfiles = 96
samples = c("BalAcu02","BalMus01","ENPAK19","ENPAK20","ENPAK21","ENPAK22","ENPAK23","ENPAK24","ENPAK25","ENPAK26","ENPAK27","ENPAK28","ENPAK29","ENPAK30","ENPBC16","ENPBC17","ENPBC18","ENPCA01","ENPCA02","ENPCA03","ENPCA04","ENPCA05","ENPCA06","ENPCA07","ENPCA08","ENPCA09","ENPOR10","ENPOR11","ENPOR12","ENPOR13","ENPWA14","ENPWA15","EubGla01","GOC002","GOC006","GOC010","GOC025","GOC038","GOC050","GOC053","GOC063","GOC068","GOC071","GOC077","GOC080","GOC082","GOC086","GOC091","GOC100","GOC111","GOC112","GOC116","GOC125","MegNov01")

sumdir = "./summary_stats/"
dir.create(sumdir, recursive = T)

# main -------
# load genotype depth files
sumdp_dtlist = lapply(1:nfiles, function(ii) {
	idx = stringr::str_pad(ii, width = 2, side = "left", pad = "0")
	myfile = paste0(dataset, "_", ref, "_chr", idx, "_", datadate, "_sumDP.tsv")
	data = read.table(file = myfile, header = TRUE, sep = "\t")
	if (!all(colnames(data)[1:length(samples)] == samples)) {
		stop("sample missing!")
	}
	return(data)
})
sumdp_dt = dplyr::bind_rows(sumdp_dtlist)

# load genotype count files
countgt_dtlist = lapply(1:nfiles, function(ii) {
	idx = stringr::str_pad(ii, width = 2, side = "left", pad = "0")
	myfile = paste0(dataset, "_", ref, "_chr", idx, "_", datadate, "_countgt.tsv")
	data = read.table(file = myfile, header = TRUE, sep = "\t")
	if (!all(colnames(data)[1:length(samples)] == samples)) {
		stop("sample missing!")
	}
	return(data)
})
countgt_dt = dplyr::bind_rows(countgt_dtlist)

# get the average for all
sumdp = sumdp_dt %>% dplyr::select(-LISTID) %>% base::colSums()
countgt = countgt_dt %>% dplyr::select(-LISTID) %>% base::colSums()
meandp = sumdp / countgt
upmeandp2 = round(2 * meandp)
upmeandp2.5 = round(2.5 * meandp)

# write the mean dp
meandp = namednum2df(meandp, c("sample", "mean_gtDP"))
write.table(meandp, file = paste0(sumdir, "all_summary_mean_gtDP.csv"),sep = ",", row.names = F, col.names = T)
# write the upmean dp
upmeandp2 = namednum2df(upmeandp2, c("sample", "upmean_gtDP"))
write.table(upmeandp2, file = paste0(sumdir, "all_summary_2xmean_gtDP.csv"),sep = ",", row.names = F, col.names = F)
upmeandp2.5 = namednum2df(upmeandp2.5, c("sample", "upmean_gtDP"))
write.table(upmeandp2.5, file = paste0(sumdir, "all_summary_2.5xmean_gtDP.csv"),sep = ",", row.names = F, col.names = F)
date()
closeAllConnections()
