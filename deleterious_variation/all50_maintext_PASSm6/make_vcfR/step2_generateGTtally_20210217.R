# Title: Generate summary of genotype tally
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Feb 18 14:52:02 2021

# preparation --------
rm(list = ls())
cat("\014")

# library(vcfR)
library(stringr)
library(dplyr)
library(ggplot2)

source("~/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# def functions --------
generate_summary <- function(gt, muttype) {
    # start with output df with sample name at each row and tally info at each column
    # CalledCount = "HomRefCount"+"HomAltCount"+"HetCount"
    # PassCount = "CalledCount"+"MissingCount"
    tallyinfo = c("SampleId","HomRefCount","HomAltCount","HetCount","CalledCount","MissingCount","PassCount")
    outdf = data.frame(matrix(nrow = ncol(gt), ncol = length(tallyinfo)))
    colnames(outdf) = tallyinfo
    outdf[,'SampleId'] = colnames(gt)
    # loop through each sample
    for(ii in 1:ncol(gt)) {
        temp = gt[,ii]
        # fix the na
        temp[is.na(temp)] = './.'
        # check if right sample
        if (outdf[ii,'SampleId']!=colnames(gt)[ii]) {
            stop('ERROR: SampleId mismatch!')
        }
        outdf[ii,'HomRefCount']=sum(temp == '0/0')
        outdf[ii, 'HomAltCount']=sum(temp == '1/1')
        outdf[ii, 'HetCount']=sum(temp == '0/1')
        outdf[ii, 'CalledCount']=sum(temp != './.')
        outdf[ii, 'MissingCount']=sum(temp == './.')
        # sanity check
        if (outdf[ii, 'CalledCount'] != outdf[ii,'HomRefCount'] + outdf[ii, 'HomAltCount'] + outdf[ii, 'HetCount']) {
            stop('ERROR: CalledCount mismatch!')
        }
        # sanity check
        if (nrow(gt) != outdf[ii,'CalledCount'] + outdf[ii, 'MissingCount']) {
            stop('ERROR: PassCount mismatch!')
        }
    }
    outdf[,'PassCount'] = nrow(gt)
    outdf = outdf %>%
        dplyr::mutate(HomRefPer = HomRefCount/CalledCount,
                      HomAltPer = HomAltCount/CalledCount,
                      HetPer = HetCount/CalledCount,
                      RefAlleleCount = 2*HomRefCount + HetCount,
                      AltAlleleCount = 2*HomAltCount + HetCount,
                      CalledAlleleCount = 2*CalledCount,
                      RefAllelePer = RefAlleleCount/CalledAlleleCount,
                      AltAllelePer = AltAlleleCount/CalledAlleleCount,
                      PopId = substr(SampleId, 1, 3),
                      SubPopId = ifelse(PopId == "ENP",
                                        substr(SampleId, 4, 5), PopId),
                      MutType = muttype)
    # sanity check: CalledAlleleCount = RefAlleleCount + AltAlleleCount
    if (!(all(outdf$CalledAlleleCount == outdf$RefAlleleCount + outdf$AltAlleleCount))) {
        stop('ERROR: Allele Count mismatch!')
    }
    return(outdf)
}

plot_summary <- function(mytally, mytitle) {
    forplot = mytally %>%
        select(ends_with('Id'), HomRefPer,HetPer,HomAltPer, MutType) %>%
        reshape2::melt(id.vars = c('SampleId', 'PopId', 'SubPopId', 'MutType'))
    forplot$variable = factor(forplot$variable, 
                              levels = c('HomRefPer','HetPer','HomAltPer'),
                              labels = c('Homozygous Ancestral', 'Heterozygous', 'Homozygous Derived'))
    pp <- ggplot(forplot, aes(x = SampleId, y = value, fill = variable)) + 
        geom_bar(position="stack", stat="identity") +
        scale_fill_brewer(palette = "Greys") +
        geom_vline(xintercept = 27.5, color = "#377EB8") +
        labs(y = "Proportion", fill = "Variant Types", title = mytitle) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
    return(pp)
}

plot_popsummary <- function(alltally) {
    forplot = alltally %>%
        select(ends_with('Id'), HomRefPer,HetPer,HomAltPer, MutType) %>%
        reshape2::melt(id.vars = c('SampleId', 'PopId', 'SubPopId', 'MutType'))
    forplot$variable = factor(forplot$variable, 
                              levels = c('HomRefPer','HetPer','HomAltPer'),
                              labels = c('Homozygous Ancestral', 'Heterozygous', 'Homozygous Derived'))
    forplot2 = forplot %>%
        dplyr::mutate(Id = paste(MutType, PopId, sep = '|')) %>%
        group_by(Id, PopId, MutType, variable) %>%
        summarise(meanval = mean(value))
    
    pp1 <- ggplot(forplot2, aes(x = Id, y = meanval, fill = variable)) + 
        geom_bar(position="stack", stat="identity") +
        scale_fill_brewer(palette = "Greys") +
        labs(y = "Proportion", fill = "Variant Types") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
    
    pp2 <- ggplot(forplot2, aes(x = MutType, y = meanval, fill = variable)) + 
        geom_bar(position="dodge", stat="identity") +
        facet_wrap(~ PopId) + 
        scale_fill_brewer(palette = "Greys") +
        labs(y = "Proportion", fill = "Variant Types") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
    
    return(list(pp1, pp2))
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
cdstype = 'ALLregions'
prefixlist = c("syn", "nonsyn", "nonsynDEL", "nonsynTOL", "LOF", "BenignRM", "Damaging")

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR/', dataset, ref, sep = '/')
setwd(workdir)

indir = './derive_data/gttable/'
plotdir = './plots/'
outdir = './derive_data/sum_table/'
dir.create(plotdir)
dir.create(outdir)

# main --------
gttallylist = lapply(prefixlist, function(xx) {
    gtseg = readRDS(file = paste0(indir, xx,'_all50_Minke_ALLregions_GTtable_Seg_PASSm6_20210218.rds'))
    mytally = generate_summary(gt = gtseg, muttype = xx)
    # output tally file 
    write.csv(mytally, file = paste0(outdir, xx, '_all50_Minke_ALLregions_SUMtable_Seg_PASSm6_', today,'.csv'))
    pp = plot_summary(mytally, mytitle = paste0(xx, '_all50_Minke_ALLregions_Seg_PASSm6'))
    ggsave(filename = paste0(xx, '_all50_Minke_ALLregions_Seg_PASSm6_', today, '.pdf'), 
           plot = pp, path = plotdir, height = 6, width = 9)
    return(mytally)
})

# summarize info ========
gttallydt = dplyr::bind_rows(gttallylist)
write.csv(gttallydt, file = paste0(outdir, 'all50_Minke_ALLregions_SUMtable_Seg_PASSm6_', today,'.csv'))
pplist = plot_popsummary(gttallydt)
ggsave(filename = paste0('all50_Minke_ALLregions_Seg_PASSm6_Summary1_', today, '.pdf'), 
       plot = pplist[[1]], path = plotdir, height = 6, width = 8)
ggsave(filename = paste0('all50_Minke_ALLregions_Seg_PASSm6_Summary2_', today, '.pdf'), 
       plot = pplist[[2]], path = plotdir, height = 6, width = 8)

# output files --------
saveRDS(gttallydt, file = paste0(outdir, 'all50_Minke_ALLregions_SUMtable_Seg_PASSm6_', today,'.rds'))

# cleanup --------
date()
closeAllConnections()
