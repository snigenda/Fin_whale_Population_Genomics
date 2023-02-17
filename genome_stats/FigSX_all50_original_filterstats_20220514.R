# Title: Plot the supplemental figure for filter statisitics
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat May 14 23:34:37 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(dplyr)
library(reshape2)
library(ggplot2)

setwd('/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/revisions_VariantFilter/previous_filterstats/')

source("~/Lab/fin_whale/scripts/config/plotting_config.R")

# def functions --------
# apply grep to multiple things
multi_grep <- function(patternlist, xx) {
    grepres <- lapply(patternlist, function(pp) {
        grep(pp,xx)
    })
    outres = unique( unlist(grepres) )
    return(outres)
}

# summarize site filter tally
summarize_sitetally <- function(df) {
    # first qual filters
    qualid = grep('FAIL_qual',df$SITE_FILTER)
    # then hard filters
    hardid = setdiff(multi_grep(c('FAIL_QD','FAIL_FS','FAIL_MQ','FAIL_MQRS', 'FAIL_RPRS', 'FAIL_SOR'),
                                df$SITE_FILTER), qualid)
    # then cpg rep
    cpgid = setdiff(grep('FAIL_CpGRep', df$SITE_FILTER), c(qualid, hardid))
    # then variant type is not SNPs
    snpid = setdiff(multi_grep(c('FAIL_mutType', 'FAIL_badRef', 'FAIL_noADi', 'FAIL_badAlt'), df$SITE_FILTER),
                    c(qualid, hardid,cpgid))
    # then excessHet
    hetid = setdiff(grep('WARN_excessHet',df$SITE_FILTER),
                    c(qualid, hardid,cpgid,snpid))
    # then missing
    missid = setdiff(multi_grep(c('FAIL_noGT', 'WARN_missing'), df$SITE_FILTER),
                     c(qualid, hardid,cpgid,snpid,hetid))
    passid = grep('PASS', df$SITE_FILTER)
    # all errors
    all_types = list(qualid, hardid,cpgid,snpid,hetid,missid,passid)
    names(all_types) = c('QUAL < 30', 'Hard Filter', 'CpG or Repeat', 'Not SNP', 'GT Het > 75%', 'Missing > 20%', 'PASS')
    
    # subset
    sumfilter = lapply(all_types,function(xx) {
        out = data.frame(t(df[xx,-1] %>% colSums()), stringsAsFactors = FALSE)
    })

    sumfilterdf = dplyr::bind_rows(sumfilter,.id = 'SITE_FILTER')
    
    # sanity check
    colSums(sumfilterdf[,-1])
    sumfilterdf$SITE_FILTER = factor(sumfilterdf$SITE_FILTER, levels = sumfilterdf$SITE_FILTER)
    return(sumfilterdf)
}

# sort by the total order 
summarize_gttally <- function(df) {
    # df keep 
    # check the rowSum should be the same as the total sites
    totalgt = colSums(apply(df[,2:51], 2, as.numeric))
    unique(totalgt) # 2324429748 (the maximum genotype called available)
    totalsites = unique(totalgt)
    
    # bin filters
    minDPid = grep('minDP',df$GT_FILTER)
    maxDPid = setdiff(grep('maxDP',df$GT_FILTER),minDPid)
    gqid = setdiff(grep('GQ20',df$GT_FILTER),c(minDPid,maxDPid))
    
    # all errors (note this DEPEND ON THE ORDERS)
    all_types = list(12,10,minDPid,maxDPid,gqid,
                     4,3,c(1,2),11)
    # sort(unlist(all_types))
    names(all_types) = c('Site FAIL', './.', 'GT_DP < 8', 'GT_DP > 2.5x mean', 'GT_GQ < 20', 
                         'GT_ABHomRef < 0.9', 'GT_ABHomAlt > 0.1', 'GT_ABHet > 0.8 or < 0.2', 'PASS')
    
    # subset
    sumfilter = lapply(all_types,function(xx) {
        out = data.frame(t(df[xx,-1] %>% colSums()), stringsAsFactors = FALSE) 
        out = out/totalsites
    })
    
    sumfilterdf = dplyr::bind_rows(sumfilter,.id = 'GT_FILTER')
    
    # sanity check
    colSums(sumfilterdf[,-1])
    sumfilterdf$GT_FILTER = factor(sumfilterdf$GT_FILTER, levels = sumfilterdf$GT_FILTER)
    return(sumfilterdf)
}


# def variables --------
# args = commandArgs(trailingOnly=TRUE)
# dataset = as.character(args[1])
# ref = as.character(args[2])
dataset = "all50"
ref = "Minke"
ncontig = 96
today = format(Sys.Date(), "%Y%m%d")
datadate = "20201110"

# load data --------
gtfiltertally = read.csv(file = '/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/Summary_stats/all50/summary_stats/gtfilter_tally_all_20201112.csv') 
gtfiltertally = gtfiltertally[,-1]
filtertally = read.csv(file = '/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/Summary_stats/all50/summary_stats/filter_tally_per_20201112.csv')
colnames(filtertally)[1] = 'SITE_FILTER'

# plot site filter --------
sitefiltersumdf = summarize_sitetally(df = filtertally) %>% 
    dplyr:: relocate(N_FILTER_all, .after = SITE_FILTER)
# add a size specifications
coorddf = data.frame(variable = colnames(sitefiltersumdf)[-1],
                              xcoord = c(0.5,5:100))
forplot1 = sitefiltersumdf %>%
    reshape2::melt(id.vars = "SITE_FILTER") %>% 
    dplyr::mutate(facet = ifelse(variable == 'N_FILTER_all', 'All', 'By contig groups')) %>%
    dplyr::left_join(., y = coorddf, by = 'variable')

pp1 <- ggplot(data = forplot1, aes(x = variable, y = value, fill = SITE_FILTER)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Contig list", y = "Site Filter %", fill = 'Site\nFilter') +
    facet_grid(. ~ facet, scales = "free_x", space = "free_x") +
    scale_x_discrete(breaks = c('N_FILTER_all', paste0('N_FILTER_', 
                                                       stringr::str_pad(seq(0,96,by=5),width=2,pad='0'))),
                     labels = c('all',seq(0,95,by=5))) + # c("all", as.character(1:96))
    scale_fill_brewer(palette = "Set3") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text = element_text(size = 9),
          axis.text.x = element_text(angle = 90, hjust = 1))

# plot genotype filters --------
gtfiltersumdf = summarize_gttally(df = gtfiltertally)
forplot2 = gtfiltersumdf %>% 
    reshape2::melt(id.vars = "GT_FILTER") 

pp2 <- ggplot(data = forplot2, aes(x = variable, y = value, fill = GT_FILTER)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "SampleId", y = "Genotype Filter %", fill = "Genotype\nFilter") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 9, name = 'Set3')[c(9,1:6,8,7)]) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() + 
    theme(axis.text = element_text(size = 9),
          axis.text.x = element_text(angle = 90, hjust = 1))

# output files --------
outpp <- ggpubr::ggarrange(pp1,pp2, nrow = 2, labels = c('A','B'))

ggsave(filename = paste0("FigSX_all50_original_filterstats", today, ".pdf"), plot = outpp, width = 11, height = 8)

write.csv(x = sitefiltersumdf, file = 'FigSX_all50_original_filterstats_site.csv')
write.csv(x = gtfiltersumdf, file = 'FigSX_all50_original_filterstats_gt.csv')

save.image(file = paste0("FigSX_all50_original_filterstats", today, ".RData"))

closeAllConnections()

# get all the number of sites from sitefiltersummary --------
# Mon Aug  8 16:22:54 2022
sitefiltercount = data.frame(2324429748 * sitefiltersumdf[,'N_FILTER_all'])
rownames(sitefiltercount) = sitefiltersumdf$SITE_FILTER
colnames(sitefiltercount) = 'N_FILTER_all'
# this might not be the most right way of calculating things. but it gives a source
sitefiltercount['PASS','N_FILTER_all'] + sitefiltercount['Missing > 20%','N_FILTER_all']
write.csv(sitefiltercount, file = 'FigSX_all50_original_filterstats_site_count.csv')



