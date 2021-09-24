# Title: generate reference genome contiglist --------
# Author: Meixi Lin
# Date: Thu Aug  6 15:35:56 2020
# Author:
# Date:
# Modification:

# preparation --------
options(echo = FALSE)
rm(list = ls())
cat("\014")

library(dplyr)
date() # the execution date

# def functions --------

# set value -------
today = format(Sys.Date(), "%Y%m%d")
genomename="GCA_009873245.2_mBalMus1.v2"
data.dir=paste0("/u/project/rwayne/snigenda/finwhale/cetacean_genomes/blue_whale_genome/", genomename)
dictname=paste0(genomename, "_genomic.dict")
assemblyname=paste0(genomename,"_assembly_report.txt")

setwd(data.dir)

# load data -------
dict <- read.delim(file = dictname, header = F, skip = 1, stringsAsFactors = F) %>%
    dplyr::tbl_df() %>% 
    dplyr::select(V2, V3, V4) %>%
    dplyr::mutate(SN = stringr::str_sub(V2, start = 4),
                  LN = as.integer(stringr::str_sub(V3, start = 4)))

assembly <- read.delim(file = assemblyname, comment.char = "", skip = 32, stringsAsFactors = F)
colnames(assembly)[1] = "Sequence.Name" 
chromosomes = assembly[assembly$Assigned.Molecule.Location.Type == "Chromosome", "GenBank.Accn"]

# # decide cutoff --------
LN <- dict$LN
sum(LN) # total genomelen
#  2380012384
table(LN > 3e+6) 
table(LN > 3e+6) 
# FALSE  TRUE 
#   108    22 
sum(LN[LN > 3e+6])/sum(LN) 
# [1] 0.9968733

# cut off --------
newdict <- dict %>% dplyr::mutate(binid = NA)

# bin the dictionary --------
binid=1
for (ii in 1:nrow(newdict)) {
# if in placed chromosome, put it in a single contiglist
    if (newdict[ii, "SN"] %in% chromosomes) {
        newdict[ii, "binid"]<-binid
        binid=binid+1
    } else {
        newdict[ii, "binid"]<-binid
    }
}

# size bin --------
sizebin <- newdict %>%
    dplyr::group_by(binid) %>%
    dplyr::summarise(sumlen = sum(LN), 
                     count = n())

# output dictionary ---------
outdict <- newdict %>% 
    dplyr::left_join(., sizebin, by = "binid") %>%
    dplyr::mutate(binid = stringr::str_pad(binid, width = 2, pad = "0")) 

# output the contig list --------
outpath = "./contiglist/"
dir.create(path = outpath, recursive = T)
write.csv(x = outdict, file = paste0(genomename, "_contig_summary_", today, ".csv"))
outlist <- dplyr::group_split(outdict, binid) 
lapply(outlist, function(xx) {
    temp <- xx$SN 
    write.table(temp, file = paste0(outpath, genomename, "_genomic.contiglist_", xx$binid[1], ".list"), quote = F, sep = "\t", row.names = F, col.names = F)
    return(0)
})
