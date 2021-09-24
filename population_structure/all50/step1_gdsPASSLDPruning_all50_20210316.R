# Title: Filter GDS file with 1. filter == PASS and 2. SNPRelated LD pruning
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin (meixilin@ucla.edu)
# Date: Wed Mar 17 11:27:31 2021
# Example: Rscript --vanilla step1_gdsLDPruning_all50_20210316.R '/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/all50/Minke' 'JointCalls_all50_08_B_VariantFiltration_bialleic_all.gds' 'JointCalls_all50_filterpass_bialleic_all'
# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(gdsfmt)
library(SNPRelate)

# def functions --------
# # parse mafcutoff
# parse_maf <- function(mafcutoff) {
#     if (mafcutoff == 'NaN') {
#         out = NaN
#     } else {
#         out = as.numeric(mafcutoff)
#     }
#     return(out)
# }

# subset by FILTER=='PASS'
subset_filter_pass <- function(genofile, gdsfile, chr.id, position.id, outprefix) {
    # CHECKED: works the same as 'bcftools view -i'TYPE="snp" & FILTER="PASS"''
    g <- (read.gdsn(index.gdsn(genofile, "snp.annot/filter")) == "PASS")
    snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))[g]
    print(length(snpset))
    # sanity check
    if(any(duplicated(snpset))) {
        stop('SNP Id not unique, check your input!')
    }

    # Write the filterpass sites
    write.table(cbind(chr.id[snpset],position.id[snpset]), file = paste0(outprefix, '_', today,'.txt'), sep = '\t', quote = FALSE, col.names = F, row.names = F)
    # Save the gds file
    outfile = paste0(outprefix,'.gds')
    snpgdsCreateGenoSet(gdsfile, outfile, snp.id=snpset)

    print(paste(date(), ': Finished subsetting FILTER==PASS. Outputting to', outfile))
    return(snpset)
}

# define function for pipeline so it can run with different methods of LD pruning
# note that it intersects with the snpset that passed filters (snpset.filt)
SNPPruning <- function(genofile, ldmethod, mafcutoff, snpset.filt) {
    ### first check the mafcutoff input
    # mafcutoff = parse_maf(mafcutoff)
    ### LD-based SNP pruning (Filter of linkage equilibrium)
    set.seed(1000)
    #LD SNP pruning: here the maf might need testing. used to be maf = 0.09 (testsix and all50 dataset). now generate by input type
    #CHECKED: the remove.mono is automatically set to be true and works the same as 'bcftools view -i 'COUNT(GT="AA")=N_SAMPLES-N_MISSING || COUNT(GT="RR")=N_SAMPLES-N_MISSING'', i.e. the missing genotypes are taken into considerations
    #EDIT: here we pass the snp.id first as the filtered snpset before performing LD pruning analyses (mattered a lot for downstream results)
    snpset <- snpgdsLDpruning(genofile, snp.id = snpset.filt, ld.threshold = 0.2, maf = mafcutoff, autosome.only = FALSE, method = ldmethod)
    snpset.id <- unlist(unname(snpset))
    return(snpset.id)
}

# run the pruning process and output file
run_SNPPruning <- function(genofile, gdsfile, chr.id, position.id, ldmethod, mafcutoff, mafprefix, outprefix, snpset.filt) {
    #Get filtered SNPs out of LD pruning with three settings
    snpset.id <-SNPPruning(genofile, ldmethod, mafcutoff, snpset.filt)

    #Write table of SNPs
    write.table(cbind(chr.id[snpset.id],position.id[snpset.id]), file = paste0(outprefix,'_LDPruned_', mafprefix, '_', today,'.txt'), sep = "\t", quote = FALSE, col.names = F, row.names = F)

    #Create new GDS file
    outfile = paste0(outprefix,'_LDPruned_', mafprefix, '.gds')
    snpgdsCreateGenoSet(gdsfile, outfile, snp.id=snpset.id)

    outmsg = paste(date(), ': Finished LDpruning using', ldmethod, 'method, with minor allele frequency cutoff maf =', mafcutoff, '. Outputting to', outfile)
    return(outmsg)
}


# def variables --------
args <- commandArgs(trailingOnly=TRUE)
workdir <- as.character(args[1]) # the working directory (should be the same as before and after)
gdsfile <- as.character(args[2]) # the name of gdsfile to analyze
outprefix <- as.character(args[3]) # the prefix of gdsfile to output
# ldmethod <- as.character(args[4]) # method to use for ld pruning
# mafcutoff <- as.character(args[5]) # maf cutoff threshold to use for ld pruning

# define ld pruning method
ldmethod = 'composite'
# define preset mafcutoff
maflist = c(NaN, 0.05, 0.1)
mafout = c('mafNA', 'maf05', 'maf10')

today = format(Sys.Date(), "%Y%m%d")

setwd(workdir)
sessionInfo()

# load data --------
snpgdsSummary(gdsfile)
# Read the gds file
genofile <- snpgdsOpen(gdsfile, readonly = TRUE)

# main --------
#chromosome id
chr.id <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
#position id: snp.position, the base position of each SNP on the chromosome, and 0 for unknown position; it does not allow NA.
position.id <- read.gdsn(index.gdsn(genofile, "snp.position"))

# Remove the filtered sites and output the filter passed files for future use
snpset.id.pass <- subset_filter_pass(genofile, gdsfile, chr.id, position.id, outprefix)

# run three set LDpruning with different maf cutoff
run_SNPPruning(genofile = genofile, gdsfile = gdsfile, chr.id = chr.id, position.id = position.id, ldmethod = ldmethod, mafcutoff = maflist[1], mafprefix = mafout[1], outprefix = outprefix, snpset.filt = snpset.id.pass)
run_SNPPruning(genofile = genofile, gdsfile = gdsfile, chr.id = chr.id, position.id = position.id, ldmethod = ldmethod, mafcutoff = maflist[2], mafprefix = mafout[2], outprefix = outprefix, snpset.filt = snpset.id.pass)
run_SNPPruning(genofile = genofile, gdsfile = gdsfile, chr.id = chr.id, position.id = position.id, ldmethod = ldmethod, mafcutoff = maflist[3], mafprefix = mafout[3], outprefix = outprefix, snpset.filt = snpset.id.pass)

# output files --------

# cleanup --------
# Close the GDS file
closefn.gds(genofile)
date()
closeAllConnections()

