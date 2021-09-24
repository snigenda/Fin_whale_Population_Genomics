# Title: check the BalAcu contig distribution --------
# Author: Meixi Lin
# Date: Thu Dec  5 19:29:57 2019
# Author:
# Date:
# Modification:

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/UCLA/Lab/fin_whale")
library(dplyr)
date() # the execution date

# set value -------
cutoff = 1e+6 
nbins = 200
indilen = 2e+7
maxlen = 3e+7

# load data -------
dict <- read.delim(file = "/Users/linmeixi/google_drive/finwhale/preprocessing/contig_list/Minke/GCF_000493695.1_BalAcu1.0_genomic.dict", header = F, skip = 1, stringsAsFactors = F) %>%
    dplyr::tbl_df() %>% 
    dplyr::select(V2, V3, V4) %>%
    dplyr::mutate(SN = stringr::str_sub(V2, start = 4),
                  LN = as.integer(stringr::str_sub(V3, start = 4)))

# # decide cutoff --------
# LN <- dict$LN
# sum(LN) # total genomelen
# total genome len 2431687698
# table(LN > 3e+6) # 193 TRUE 
# sum(LN[LN > 3e+6])/sum(LN) # 88.7%

# table(LN > 1e+6) # 284 TRUE
# sum(LN[LN > 1e+6])/sum(LN) # 95.6% 

# table(LN > 2e+7) # 27 TRUE 
# sum(LN[LN > 2e+7])/sum(LN) # 31.4%

# cut off --------
newdict <- dict %>% 
    dplyr::filter(LN > 1e+6) %>%
    dplyr::mutate(binid = NA)
newLN <- newdict$LN
hist(newLN)

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
plot(sizebin$sumlen)
abline(h = indilen, col = "green")

# get the indibin ---------
indibin <- newdict %>%
    dplyr::filter(LN >= indilen)
idbin <- sizebin %>% 
    dplyr::filter(binid %in% indibin[['binid']])

# output dictionary ---------
outdict <- newdict %>% 
    dplyr::left_join(., sizebin, by = "binid") %>%
    dplyr::mutate(binid = stringr::str_pad(binid, width = 2, pad = "0")) 

# output the contig list --------
outpath = "./contig_list/Minke/contiglist/"
dir.create(path = outpath, recursive = T)
write.csv(x = outdict, file = paste0(outpath, "minke_contig_summary.csv"))
outlist <- dplyr::group_split(outdict, binid) 
lapply(outlist, function(xx) {
    temp <- xx$SN 
    write.table(temp, file = paste0(outpath, "BalAcu1.0_genomic.contiglist_", xx$binid[1], ".list"), quote = F, sep = "\t", row.names = F, col.names = F)
    return(0)
})
