# Title: Convert vcf files to gds files
# Author: Meixi Lin (meixilin@ucla.edu), modified by Paulina Nunez Valencia (pnunez@lcg.unam.mx)
# Date: Tue Feb 23 20:39:40 2021
# Example: Rscript --vanilla step0_vcf2gds_f50b4_20210223.R '/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/filteredvcf/f50b4/Minke' 'JointCalls_f50b4_08_B_VariantFiltration' '/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke'

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(gdsfmt)
library(SNPRelate)

# def functions --------

# def variables --------
args <- commandArgs(trailingOnly=TRUE)
vcfdir <- as.character(args[1]) # the directory of the vcffile
vcfprefix <- as.character(args[2]) # pattern for files
outdir <- as.character(args[3]) # where to output files

# vcfdir='/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/filteredvcf/f50b4/Minke'
# vcfprefix='JointCalls_f50b4_08_B_VariantFiltration'
# outdir='/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke'

dir.create(outdir, recursive = TRUE)
setwd(outdir)
sessionInfo()

# load data --------
#Get all the vcf files
vcffiles <- list.files(path = vcfdir, pattern= vcfprefix, full.names = TRUE)
vcffiles <- grep(".gz$", vcffiles, value = TRUE) # only keep the vcf.gz files
out.fn <- paste0(vcfprefix, "_bialleic_all.gds")

print(vcffiles)

# get the population names
popmap = read.csv(file = '/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/baleen_popid.csv', stringsAsFactor = FALSE)

# main --------
# convert the file
# snp.rs.id, a character string for reference SNP ID that may not be unique.
# snp.allele, it is not necessary for the analysis, but it is necessary when merging genotypes from different platforms. The format of snp.allele is “A allele/B allele”, like “T/G” where T is A allele and G is B allele.
# There are four possible values stored in the variable genotype: 0, 1, 2 and 3. For bi-allelic SNP sites, “0” indicates two B alleles, “1” indicates one A allele and one B allele, “2” indicates two A alleles, and “3” is a missing genotype. Note: "0" --> "1/1"; "1" --> "0/1"; "2" --> "0/0"; "3" --> "./."
# The SNP matrix is nsample x nsnp (snpfirstdim=FALSE) by default. Samples are in rows and SNPs are in columns. Different from VCF file specifications.
# The SNPs being taken in is not only passing sites, but all the biallelic sites (variants). Invariant sites that passed all filtering are removed. So any sites in column 5 that is A/T/C/G is included.

snpgdsVCF2GDS(vcffiles, out.fn, method = "biallelic.only")
snpgdsSummary(out.fn)

# open this file
genofile <- snpgdsOpen(out.fn, readonly = FALSE)
head(genofile)
# get the genotype data for a brief view, here, it's the first five samples and the first three sites (CHECK)
g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(5,3));g

# get the sample id
sampleid <- read.gdsn(index.gdsn(genofile, "sample.id"))
# check if the sample id matches the popmap
if (!all(sampleid == popmap$SampleId)) {
    stop('Wrong popmap file!')
}

# add population annotation
samp.annot <- data.frame(pop.group = popmap$PopId)
add.gdsn(genofile, "sample.annot", samp.annot)

# Close the GDS file
closefn.gds(genofile)

# cleanup --------
date()
closeAllConnections()
