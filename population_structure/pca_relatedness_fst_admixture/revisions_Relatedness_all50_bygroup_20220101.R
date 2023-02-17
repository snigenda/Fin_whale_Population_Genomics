# Title: revisions for kinship within major groups
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Jan  1 16:41:41 2022
# Modification: Add source_data
# Date: Sun Jan 22 16:10:40 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(dunn.test)

source('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/plotting_config.R')

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

setwd('/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/PopStructure/all50/')
outdir = './derive_data/revisions_Kinship/'
plotdir = './plots/revisions_Kinship/'

dir.create(outdir)
dir.create(plotdir)

dataset = 'all50'
ref = 'Minke'
mafcut = '10' # mafcutoff across the population used in the input gdsfile (NOT the maf setting used in the other samples)
missing = 0.2
gdsfile = 'JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds'

sessionInfo()
getwd()

# color palettes ========
subpop2_cols = loccolors[c('AK','OR','CA','GOC')]
names(subpop2_cols)[2] = 'MENP'
# subpop2_cols[2] = RColorBrewer::brewer.pal(n = 8, 'Set1')[8]
# pals::pal.bands(subpop2_cols)

# def function --------
estimate_kinship <- function(genofile, grouplev = c("PopId","SubPopId","SubPopId2"), targetpop) {
    sampleids = popmap[popmap[,grouplev] == targetpop,'SampleId']
    # Estimating IBD Using PLINK method of moments (MoM)
    # Missing rate and MAF set as the same in the tutorial (http://corearray.sourceforge.net/tutorials/SNPRelate/#f_st-estimation) and same as Paulina's previous settings
    print(paste0('Group by = ', grouplev, '. Target population = ', targetpop))
    ibddt <- snpgdsIBDMoM(genofile, sample.id = sampleids, verbose = TRUE,
                           maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE)
    # k0	the probability of sharing ZERO alleles
    # k1	the probability of sharing ONE alleles
    # kinship	kinship coefficient
    ibddt.coeff <- snpgdsIBDSelection(ibddt)
    ibddt.coeff = dplyr::left_join(ibddt.coeff, popmap, by = c('ID1' = 'SampleId'))
    return(ibddt.coeff)
}

reorder_matrix <- function(ibd.coeff) {
    ibdnames = c('ID1','ID2','kinship')
    newibd = data.frame(matrix(nrow=nrow(ibd.coeff),ncol=length(ibdnames)))
    colnames(newibd) = ibdnames
    rowid=1
    for (ii in 1:(nrow(popmap)-1)) {
        for (jj in (ii+1):nrow(popmap)) {
            newibd[rowid,'ID1']=popmap[ii,'SampleId']
            newibd[rowid,'ID2']=popmap[jj,'SampleId']
            rawid1 = (ibd.coeff$ID1 == newibd[rowid,'ID1'] & ibd.coeff$ID2 == newibd[rowid,'ID2'])
            rawid2 = (ibd.coeff$ID2 == newibd[rowid,'ID1'] & ibd.coeff$ID1 == newibd[rowid,'ID2'])
            if (sum(rawid1) == 1 & sum(rawid2) == 0) {
                newibd[rowid,'kinship'] = ibd.coeff[rawid1, 'kinship']
            } else {
                if (sum(rawid1) == 0 & sum(rawid2) == 1) {
                    newibd[rowid,'kinship'] = ibd.coeff[rawid2, 'kinship']
                }
            }
            rowid=rowid+1
        }
    }
    newibd = newibd %>%
        dplyr::filter(!is.na(kinship)) %>%
        dplyr::left_join(., popmap, by = c('ID1' = 'SampleId'))
    return(newibd)
}

plot_kinship <- function(ibd.coeff, grouplev = c("PopId","SubPopId","SubPopId2"), mycolors, reorder=FALSE) {
    # find another order if needed
    if (reorder==TRUE) {
        ibd.coeff <- reorder_matrix(ibd.coeff)
    }
    # prepare data
    ibd.coeff$ID1 = factor(ibd.coeff$ID1, levels = popmap$SampleId)
    ibd.coeff$ID2 = factor(ibd.coeff$ID2, levels = popmap$SampleId)
    ibd.coeff[,grouplev] = factor(ibd.coeff[,grouplev], levels = unique(popmap[,grouplev]))
    pp1 <- ggplot(ibd.coeff, aes(x=ID1, y=ID2, fill=kinship)) +
        geom_tile() + # aes(color = !!sym(grouplev))
        theme_light() +
        scale_fill_gradient(low = '#1F618D', high = '#2ECC71') +
        # scale_color_manual(values = mycolors) +
        theme(axis.text.x=element_text(angle=90),
              aspect.ratio = 1)
    # get a boxplot
    pp2 <- ggplot(ibd.coeff, aes(x = !!sym(grouplev), y = kinship, color = !!sym(grouplev))) +
        geom_boxplot() +
        scale_color_manual(values = mycolors) +
        # ggpubr::stat_compare_means(method = 'kruskal.test', size = 3, label.x = 1) +
        theme_light() +
        theme(legend.position = 'none',
              axis.title.x = element_blank(),
              aspect.ratio = 1)
    # set for inset
    pp <- pp1 +
        annotation_custom(
            ggplotGrob(pp2),
            xmin = 24, xmax = 44, ymin = 1, ymax = 21
        )
    return(pp)
}

get_kincompare <- function(dt1, dt2, dt3) {
    datalist <- lapply(list(dt1, dt2, dt3), function(dt) {
        data = tibble(dt[,c('ID1','ID2','kinship')])
    })
    # senp: separated enp; menp: merged enp (BC+OR+WA)
    outdt = dplyr::full_join(datalist[[1]],datalist[[2]],by=c('ID1','ID2'), suffix = c('.senp','.menp'))
    outdt = dplyr::full_join(outdt,datalist[[3]],by=c('ID1','ID2')) %>%
        dplyr::rename(kinship.enp = kinship) %>%
        dplyr::mutate(enp.senp = kinship.enp - kinship.senp,
                      enp.menp = kinship.enp - kinship.menp)
    # total should be 30*29/2 + 20*19/2 = 625 rows
    # join popmap
    outdt = dplyr::left_join(outdt, popmap, by = c('ID1' = 'SampleId'))
    return(outdt)
}

kinship_stats <- function(ibd, grouplev = c("PopId","SubPopId","SubPopId2")) {
    print('#######################################################################')
    # excluding zeros
    ibdsum = ibd %>%
        dplyr::group_by(!!sym(grouplev)) %>%
        dplyr::summarise(mean_kinship = mean(kinship),
                         mean_kinship_plus0 = sum(kinship)/sum(kinship>0),
                         .groups = 'drop')
    print(ibdsum)
    print(ibd[,c(grouplev,'kinship')] %>% split(., ibd[,grouplev]) %>% purrr::map(summary))

    # run dunn tests and kruskal tests
    ibd[,grouplev] = factor(ibd[,grouplev], levels = unique(popmap[,grouplev]))
    print(kruskal.test(ibd[,'kinship'] ~ ibd[,grouplev]))
    if (nlevels(ibd[,grouplev]) > 2) {
        dunn.test::dunn.test(x = ibd[,'kinship'], g = ibd[,grouplev], method = 'none')
        dunn.test::dunn.test(x = ibd[,'kinship'], g = ibd[,grouplev], method = 'bonferroni')
        dunn.test::dunn.test(x = ibd[,'kinship'], g = ibd[,grouplev], method = 'bh')
    }
    return(ibdsum)
}

# load data --------
genofile = SNPRelate::snpgdsOpen(filename = paste0(ref, '/', gdsfile))

#List of sample ids
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# sample designations
popmap = read.csv(file = '/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/popmap_all50.csv', stringsAsFactors = FALSE) %>%
    dplyr::mutate(SubPopId2 = ifelse(SubPopId %in% c('WA', 'OR', 'BC'), 'MENP', SubPopId)) %>%
    dplyr::mutate(SubPopId2 = factor(SubPopId2, levels = c('AK','MENP','CA','GOC'))) %>%
    dplyr::arrange(SubPopId2)

# main --------
# estimate IBD within sub population ========
subpop_ibdlist <- lapply(mylocs, function(xx) {
    ibd <- estimate_kinship(genofile, grouplev = "SubPopId", targetpop = xx)
})

subpop_ibd <- dplyr::bind_rows(subpop_ibdlist)

pp1 <- plot_kinship(subpop_ibd, mycolors = loccolors, grouplev = "SubPopId")
ggsave(filename = paste0(plotdir, 'kinship_subpop_', today, '.pdf'), plot = pp1, height = 8, width = 8)

# get summary ########
subpop_mean <- kinship_stats(ibd = subpop_ibd, grouplev = 'SubPopId')

# estimate IBD within MENP groups ========
subpop2_ibdlist <- lapply(c('AK','CA','MENP','GOC'), function(xx) {
    ibd <- estimate_kinship(genofile, grouplev = "SubPopId2", targetpop = xx)
})

subpop2_ibd <- dplyr::bind_rows(subpop2_ibdlist)

pp2 <- plot_kinship(subpop2_ibd, mycolors = subpop2_cols, grouplev = "SubPopId2")
ggsave(filename = paste0(plotdir, 'kinship_subpop2_', today, '.pdf'), plot = pp2, height = 8, width = 8)

# get summary ########
subpop2_mean <- kinship_stats(ibd = subpop2_ibd, grouplev = 'SubPopId2')

# estimate using two groups ========
pop_ibdlist <- lapply(mypops, function(xx) {
    ibd <- estimate_kinship(genofile, grouplev = "PopId", targetpop = xx)
})

pop_ibd <- dplyr::bind_rows(pop_ibdlist)

pp3 <- plot_kinship(pop_ibd, mycolors = mycolors, grouplev = "PopId", reorder = TRUE)
ggsave(filename = paste0(plotdir, 'kinship_pop_', today, '.pdf'), plot = pp3, height = 8, width = 8)

# get summary ########
pop_mean <- kinship_stats(ibd = pop_ibd, grouplev = 'PopId')

# get enp comparisons ========
allkinships <- get_kincompare(subpop_ibd, subpop2_ibd, pop_ibd)
write.csv(allkinships, file = paste0(outdir, 'ENPsub_kinships_', today, '.csv'))

# cleanup --------
closefn.gds(genofile)
date()

# save image --------
save.image(file = paste0(outdir, 'revisions_Relatedness_all50_bygroup_', today, '.RData'))
closeAllConnections()

# load the image to output the data --------
load('./derive_data/revisions_Kinship/revisions_Relatedness_all50_bygroup_20220102.RData')
# remove the comparisons
outplotdata = allkinships %>%
    dplyr::select(-enp.senp,-enp.menp)

# rename the headers
colnames(outplotdata)[3:5] = c('kinship_B_separate_ENP', 'kinship_C_middle_ENP', 'kinship_A_all_ENP')
write.csv(outplotdata,file = '~/Lab/fin_whale/FinWhale_PopGenomics_2021/source_data/FigS11.csv')
