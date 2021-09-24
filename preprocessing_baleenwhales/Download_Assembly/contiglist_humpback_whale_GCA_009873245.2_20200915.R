# Title: generate reference genome contiglist --------
# Author: Meixi Lin
# Date: "Mon Oct 12 11:31:57 2020"
# Author:
# Date:
# Modification:

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

library(dplyr)
date() # the execution date

# def functions --------

# set value -------
cutoff = 1e+6 
nbins = 200
indilen = 2e+7
maxlen = 4e+7

today = format(Sys.Date(), "%Y%m%d")
genomename="GCA_004329385.1_megNov1"
data.dir=paste0("/u/project/rwayne/snigenda/finwhale/cetacean_genomes/humpback_whale_genome/", genomename)
# # local 
# data.dir=paste0("~/google_drive/finwhale/cetacean_genomes/humpback_whale_genome/", genomename)
dictname=paste0(genomename, "_genomic.dict")
assemblyname=paste0(genomename,"_assembly_report.txt")

setwd(data.dir)

# load data -------
dict <- read.delim(file = dictname, header = F, skip = 1, stringsAsFactors = F) %>%
    tibble::as_tibble() %>% 
    dplyr::select(V2, V3, V4) %>%
    dplyr::mutate(SN = stringr::str_sub(V2, start = 4),
                  LN = as.integer(stringr::str_sub(V3, start = 4)))

assembly <- read.delim(file = assemblyname, comment.char = "", skip = 31, stringsAsFactors = F)
colnames(assembly)[1] = "Sequence.Name" 
table(assembly$Sequence.Role) # all unplaced-scaffold

# # decide cutoff --------
LN <- dict$LN
sum(LN) # total genomelen
#  2265788366
table(LN > 3e+6) 
# FALSE  TRUE 
# 2320   238 
table(LN > 1e+6) 
# FALSE  TRUE 
# 2197   361 
sum(LN[LN > 3e+6])/sum(LN) 
# [1] 0.8603335
sum(LN[LN > 1e+6])/sum(LN) 
# [1] 0.9626988

# cut off --------
newdict <- dict %>% 
    dplyr::filter(LN > 1e+6) %>%
    dplyr::mutate(binid = NA) 
newdict$binid = as.integer(newdict$binid)
sum(newdict$LN)
# hist((newdict$L)

# get all the other shorter fragments 
shortdict <- dict %>% 
    dplyr::filter(LN <= 1e+6)
sum(shortdict$LN)
# hist(shortdict$LN)

# bin the dictionary --------
for (ii in 1:nbins) {
    # select the top one not binned yet 
    if (any(is.na(newdict$binid)) == FALSE) {
        break;
    }
    index = which(is.na(newdict$binid))[1]
    temp = as.integer(newdict[index, "LN"])
    newdict[index, "binid"] = ii
    if (temp >= indilen) {
        index = index + 1
    } 
    else {
        while (temp <= maxlen & index <= nrow(newdict)) {
            testid = index + 1
            if (testid > nrow(newdict)) {
                break;
            }
            testlen = as.integer(newdict[testid, "LN"])
            if (testlen >= indilen) {
                index = index + 1
                break;
            }
            testtemp = temp + testlen
            if (testtemp > maxlen) {
                index = index + 1
                break;
            }
            index = index + 1
            temp = testtemp
            newdict[index, "binid"] = ii
        }
    }
}

# size bin --------
sizebin <- newdict %>%
    dplyr::group_by(binid) %>%
    dplyr::summarise(sumlen = sum(LN), 
                     count = n()) 
# plot(sizebin$sumlen)
# abline(h = indilen, col = "green")

# output dictionary ---------
outdict <- newdict %>% 
    dplyr::left_join(., sizebin, by = "binid") %>%
    dplyr::mutate(binid = stringr::str_pad(binid, width = 2, pad = "0")) 

# get the indibin ---------
indidict <- outdict %>%
    dplyr::filter(LN >= indilen)

# output the contig list --------
outpath = "./contiglist/"
dir.create(path = outpath, recursive = T)
write.csv(x = outdict, file = paste0(genomename, "_contig_summary_", today, ".csv"))
write.csv(x = indidict, file = paste0(genomename, "_contig_summary_L2e7_INDI_", today, ".csv"))
write.csv(x = shortdict, file = paste0(genomename, "_contig_summary_l1e6_NOTINCLUDED_", today, ".csv"))
outlist <- dplyr::group_split(outdict, binid) 
lapply(outlist, function(xx) {
    temp <- xx$SN 
    write.table(temp, file = paste0(outpath, genomename, "_genomic.contiglist_", xx$binid[1], ".list"), quote = F, sep = "\t", row.names = F, col.names = F)
})

