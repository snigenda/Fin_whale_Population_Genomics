# Title: Convert the LD pruned gds files to vcf files
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Mon Mar  1 01:29:50 2021
# Example: Rscript --vanilla  '/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA.gds' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA'

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)

# def functions --------

# def variables --------
args <- commandArgs(trailingOnly=TRUE)
workdir <- as.character(args[1]) # the working directory (should be the same as before and after)
gdsfile <- as.character(args[2]) # the name of gdsfile to analyze
outprefix <- as.character(args[3]) # the prefix of vcffile to output
# majorref <- as.logical(args[4]) # should you use the major allele as a reference allele (USUALLY SET AS FALSE)
majorref <- FALSE
today = format(Sys.Date(), "%Y%m%d")

setwd(workdir)
sessionInfo()

# load data --------

# main --------
# note this: https://github.com/zhengxwen/SeqArray/issues/22
# default to: major.ref	= TRUE
# if TRUE, use the major allele as a reference allele; otherwise, use A allele in SNP GDS file as a reference allele
outfile = paste0(outprefix,"_SeqArray_majorref", majorref,".gds")
SeqArray::seqSNP2GDS(gdsfile, outfile, verbose=TRUE, optimize=FALSE, major.ref = majorref)

# open the file
SeqArray::seqSummary(outfile)
genofile <- SeqArray::seqOpen(outfile)

# convert to vcf
SeqArray::seqGDS2VCF(genofile, vcf.fn = paste0(outprefix,"_majorref", majorref,".vcf.gz"))

# close the GDS file
seqClose(genofile)

# output files --------

# cleanup --------
date()
closeAllConnections()

