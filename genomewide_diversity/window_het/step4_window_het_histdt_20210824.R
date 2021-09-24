# Title: Generate histogram data from the `all50` and `f50b4` dataset
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Aug 24 20:31:15 2021

# Output: table histogram tally of percent of windows contained low heterozygosity values
# Relevant values: (main text) In GOC individuals ... increased average proportion of genomic regions with low heterozygosity (46% of windows contain < 1 het/kb). In contrast, the ENP population ... few regions of low heterozygosity (12% of windows with < 1 het/kb).
# Relevant values: (Figure S6 caption) The GOC individuals have a higher proportion of windows with null or very low heterozygosity (0 – 0.1 het/kb) compared with any of the other whale species, representing 13.3% – 18.6% of the total number of windows for this population.

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
# def functions --------
# raw histogram method for any given vectors (left closed)
raw_hist <- function(hetpkb, breaks) {
    hetcat = vector(length = length(breaks)-1)
    for (ii in 1:length(breaks)-1) {
        jj = ii + 1
        hetcat[ii] = sum(hetpkb >= breaks[ii] & hetpkb < breaks[jj])
    }
    return(hetcat)
}

# nwin_all50 = 2284; nwin_f50b4 = 2258
get_histdt <- function(plotdt, histbreaks, nwin) {
    sample = unique(plotdt$sample)
    nsample = length(sample)
    nbins = length(histbreaks)-1
    histdt = data.frame(matrix(ncol = nsample+2, nrow = nbins))
    colnames(histdt) = c('minhetpkb_in', 'maxhetpkb_ex', sample)
    histdt$minhetpkb_in = histbreaks[-length(histbreaks)]
    histdt$maxhetpkb_ex = histbreaks[-1]
    for (ii in sample) {
        hetpkb = plotdt[plotdt$sample == ii,'hetpkb']
        if (length(hetpkb) != nwin) {
            stop('Wrong hetpkb partition')
        }
        hetcat = raw_hist(hetpkb = hetpkb, breaks = histbreaks)
        # summary of hetcat should be nwin
        histdt[,ii] = hetcat
    }
    return(histdt)
}

# def variables --------
setwd('~/Lab/finwhale_manuscript/')

outdir = './data/window_het/derived_data/'
nwin_all50 = 2284
nwin_f50b4 = 2258

today = format(Sys.Date(), "%Y%m%d")

# load data --------
all50dt = readRDS(file = './data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_20210824.rds')

# main --------
# all50 dataset ========
# partition by 0.1
all50hist01 = get_histdt(plotdt = all50dt, 
                         histbreaks = seq(0,10,by=0.1), 
                         nwin = nwin_all50)

# get the percentage
all50per01 = cbind(all50hist01[,1:2], all50hist01[,3:ncol(all50hist01)]/nwin_all50)

# percent of windows < 1 het/kb
perlt1dt = as.data.frame(colSums(all50hist01[1:10,3:52])/nwin_all50) %>%
    tibble::rownames_to_column(var = 'sample')
colnames(perlt1dt)[2] = 'Percent_lt_1hetpkb'
# by population
mean(perlt1dt[stringr::str_detect(perlt1dt$sample,'ENP'),'Percent_lt_1hetpkb'])
# [1] 0.1199066
mean(perlt1dt[stringr::str_detect(perlt1dt$sample,'GOC'),'Percent_lt_1hetpkb'])
# [1] 0.4602671

# percent of windows < 0.1 het/kb in GOC (removing the GOC010)
gocper01 = all50per01 %>% 
    dplyr::select(starts_with('GOC')) %>%
    dplyr::select(-GOC010) %>%
    t()

min(gocper01[,1])
# [1] 0.132662
max(gocper01[,1])
# [1] 0.1860771

# output files --------
write.csv(all50hist01, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_all50_histdt_', today, '.csv'))
write.csv(all50per01, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_all50_perdt_', today, '.csv'))
write.csv(perlt1dt, file = paste0(outdir, 'winHet_1Mbwin_1Mbstep_20Per_all50_perlt1_', today, '.csv'))

# cleanup --------
date()
closeAllConnections()
