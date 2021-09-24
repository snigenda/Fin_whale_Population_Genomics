# Title: Read genotypes and calculate Rxy for all50 dataset with bad individuals removed
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Mar  5 09:40:40 2021
# Modification: Change the definition for R2xy
# Date: Fri Mar 19 10:54:20 2021


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

# library(vcfR)
library(stringr)
library(dplyr)
library(ggplot2)

# source("~/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# def functions --------
# get a tally
tally_site <- function(site, ploidy = 2) {
    # get sample n_x or n_y (haploid)
    nn = length(site[!is.na(site)]) * ploidy
    # get d_x or d_y count of mutant allele at that site
    c_het = sum(site == '0/1', na.rm = TRUE)
    c_homalt = sum(site == '1/1', na.rm = TRUE)
    dd = ploidy*c_homalt + c_het
    output = c(nn,dd)
    names(output) = c('nn', 'dd')
    return(output)
}

calc_lixy <- function(dnx, dny) {
    lixy = (dnx['dd']/dnx['nn'])*(1-dny['dd']/dny['nn'])
    lixy = unname(lixy)
    return(lixy)
}

calc_li2xy <- function(dnx, dny) {
    ratiox = (2*dnx['dd']*(dnx['nn']-dnx['dd']))/(dnx['nn']*(dnx['nn']-1))
    ratiox = unname(ratiox)
    ratioy = (2*dny['dd']*(dny['nn']-dny['dd']))/(dny['nn']*(dny['nn']-1))
    ratioy = unname(ratioy)
    # NOTE: modification 
    li2xy = (1-ratiox)*ratioy
    return(li2xy)
}

# genotype should have dimnames as the sample name
# output per site lixy or li2xy
get_lxy <- function(gtseg, xpop, ypop) {
    popx = str_subset(dimnames(gtseg)[[2]], xpop)
    popy = str_subset(dimnames(gtseg)[[2]], ypop)
    # intialize the list
    outLdt = data.frame(matrix(nrow = nrow(gtseg), ncol = 4))
    colnames(outLdt) = c('Lxnoty', 'Lynotx', 'L2xnoty', 'L2ynotx')

    # loop through all sites (use the Minke whale as the outgroup)
    for (ii in 1:nrow(gtseg)) {
        # check that this site is a variable site
        site = gtseg[ii,]
        temp = site[!is.na(site)]
        if(all(temp == '1/1') | all(temp == '0/0')) {
            stop('Not a variable site')
        }
        # get the population x and y sites
        sitex = site[popx]
        sitey = site[popy]
        # tally up the sites
        dnx = tally_site(sitex)
        dny = tally_site(sitey)
        # calculate Lixy
        li_xnoty = calc_lixy(dnx,dny)
        li_ynotx = calc_lixy(dny,dnx)
        # calculate Li2xy
        li2_xnoty = calc_li2xy(dnx,dny)
        li2_ynotx = calc_li2xy(dny,dnx)
        # append each line's results
        outLdt[ii, 'Lxnoty'] = li_xnoty
        outLdt[ii, 'Lynotx'] = li_ynotx
        outLdt[ii, 'L2xnoty'] = li2_xnoty
        outLdt[ii, 'L2ynotx'] = li2_ynotx
    }
    return(outLdt)
}

# get RRxy
calc_RRxy <- function(lxydt) {
    LLxnoty = sum(lxydt$Lxnoty)
    LLynotx = sum(lxydt$Lynotx)
    LL2xnoty = sum(lxydt$L2xnoty)
    LL2ynotx = sum(lxydt$L2ynotx)
    RRxy = LLxnoty/LLynotx
    RR2xy = LL2xnoty/LL2ynotx
    output = c(RRxy, RR2xy)
    names(output) = c('RRxy', 'RR2xy')
    return(output)
}

# FIX: I AM NOT SURE jackknife estimate
# https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/lecture2.pdf
# https://si.biostat.washington.edu/sites/default/files/modules/2017_sisg_1_9_v3_0.pdf

stderr_jackknife <- function(x) {
    n = length(x)
    correction = (n-1)/n
    varn = sum((x-mean(x))^2)
    sigmajk = sqrt(correction*varn)
    return(sigmajk)
}

# get quantile 95%
quantile_jackknife <- function(x, probs = c(0.025, 0.975)) {
    output = quantile(x, probs = probs)
    return(output)
}

# format the jackknife estimates
format_CI <- function(x) {
    # both quantile and stderr estimates
    meanx = mean(x)
    sigx = stderr_jackknife(x)
    low2sig = meanx - 2*sigx
    up2sig = meanx + 2*sigx
    outnames = c('low025', 'up975', 'low2sig', 'up2sig')
    output = c(quantile_jackknife(x), low2sig, up2sig)
    names(output) = outnames
    return(output)
}

# run block jackknife
block_jackknife <- function(fulldt, nblock = 1000) {
    blocksize = as.integer(ceiling(nrow(fulldt)/nblock))
    # leave one block out
    jksummary = data.frame(matrix(nrow = nblock, ncol = 4))
    colnames(jksummary) = c('runNum', 'blocksize', 'RRxy', 'RR2xy')
    for (ii in 1:nblock) {
        rmindex = seq((ii-1)*blocksize+1, ii*blocksize)
        jkindex = setdiff(1:nrow(fulldt), rmindex)
        # print(paste('Running Jackknife block', ii, ': block size =', length(jkindex)))
        jkdt = fulldt[jkindex,]
        jkstats = calc_RRxy(jkdt)
        jksummary[ii, ] = c(ii, length(jkindex), jkstats['RRxy'], jkstats['RR2xy'])
    }
    # calculate the sigma
    sig_RRxy = stderr_jackknife(jksummary$RRxy)
    sig_RR2xy = stderr_jackknife(jksummary$RR2xy)
    CI_RRxy = format_CI(jksummary$RRxy)
    CI_RR2xy = format_CI(jksummary$RR2xy)
    rxy_real = calc_RRxy(fulldt)
    # append the real value to jksummary (note that this has to come after sig_RR calculations)
    jksummary[nblock+1,] = c(-1, nrow(fulldt), rxy_real['RRxy'], rxy_real['RR2xy'])
    # format summary stats
    rxy_summary = c(rxy_real['RRxy'], sig_RRxy, CI_RRxy)
    r2xy_summary = c(rxy_real['RR2xy'], sig_RR2xy, CI_RR2xy)
    names(rxy_summary) = c('Rxy', 'se_Rxy','low025_Rxy', 'up975_Rxy', 'low2sig_Rxy', 'up2sig_Rxy')
    print(rxy_summary)
    names(r2xy_summary) = c('R2xy', 'se_R2xy','low025_R2xy', 'up975_R2xy', 'low2sig_R2xy', 'up2sig_R2xy')
    print(r2xy_summary)
    output = list(rxy_summary, r2xy_summary, jksummary)
    return(output)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
cdstype = 'ALLregions'
prefixlist = c("syn", "nonsyn", "nonsynDEL", "nonsynTOL", "LOF", "BenignRM", "Damaging")

xpop = 'GOC'
ypop = 'ENP'

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/Rxy', dataset, ref, sep = '/')
dir.create(workdir, recursive = TRUE)
setwd(workdir)

indir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR', dataset, ref, 'derive_data/gttable/', sep = '/')
plotdir = './plots/'
outdir1 = './derive_data/lxy_table/'
outdir2 = './derive_data/rxy_jktable/'
outdir3 = './derive_data/rxy_table/'
dir.create(plotdir)
dir.create(outdir1, recursive = T)
dir.create(outdir2, recursive = T)
dir.create(outdir3, recursive = T)
sessionInfo()

# main --------
rxyoutputlist = lapply(prefixlist, function(xx) {
    print(paste0(date(), ': Calculating Lxy on ', paste0(indir, xx,'_all50_Minke_ALLregions_GTtable_Seg_PASSm6_20210218.rds')))
    gtseg = readRDS(file = paste0(indir, xx,'_all50_Minke_ALLregions_GTtable_Seg_PASSm6_20210218.rds'))
    lxydt = get_lxy(gtseg, xpop = 'GOC', ypop = 'ENP')
    saveRDS(lxydt, file = paste0(outdir1, xx, '_all50_Minke_ALLregions_Lxy_xGOCyENP_PASSm6_', today, '.rds'))
    # get the summary and jackknife
    rxy_jk = block_jackknife(lxydt)
    # append the MutType info
    rxy_jk[[1]] = data.frame(t(rxy_jk[[1]])) %>% dplyr::mutate(MutType = xx)
    rxy_jk[[2]] = data.frame(t(rxy_jk[[2]])) %>% dplyr::mutate(MutType = xx)
    rxy_jk[[3]] = rxy_jk[[3]] %>% dplyr::mutate(MutType = xx)
    # output jackknife file individually
    saveRDS(rxy_jk[[3]], file = paste0(outdir2, xx, '_all50_Minke_ALLregions_RxyJK1000_xGOCyENP_PASSm6_', today, '.rds'))
    return(rxy_jk)
})
names(rxyoutputlist) = prefixlist

# output the summary Rxy statistics ========
rxy_summary = dplyr::bind_rows(lapply(rxyoutputlist, function(x) {x[[1]]}))
r2xy_summary = dplyr::bind_rows(lapply(rxyoutputlist, function(x) {x[[2]]}))

saveRDS(rxy_summary, file = paste0(outdir3, 'Summary_all50_Minke_ALLregions_Rxy_xGOCyENP_PASSm6_', today, '.rds'))
saveRDS(r2xy_summary, file = paste0(outdir3, 'Summary_all50_Minke_ALLregions_R2xy_xGOCyENP_PASSm6_', today, '.rds'))

# output the summary jackknife ========
jk_forplot = dplyr::bind_rows(lapply(rxyoutputlist, function(x) {x[[3]]}))
saveRDS(jk_forplot, file = paste0(outdir3, 'Summary_all50_Minke_ALLregions_1000JK_xGOCyENP_PASSm6_', today, '.rds'))

# QC plotting ========
forplot = jk_forplot %>%
    dplyr::filter(runNum > 0) %>%
    dplyr::select(MutType, RRxy, RR2xy) %>%
    reshape2::melt(id.vars = 'MutType')

forplot$variable = factor(forplot$variable, levels = c('RRxy', 'RR2xy'), labels = c('Rxy', 'R2xy'))
pp <- ggplot(data = forplot, aes(x = MutType, y = value, color = MutType)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, color = 'darkgray', linetype = 'dashed') +
    facet_wrap(. ~  variable) +
    coord_flip() +
    ggpubr::theme_pubr()

ggsave(filename = paste0('QC_Rxy_all50_Minke_ALLregions_1000JK_xGOCyENP_PASSm6.pdf'),
       plot = pp, path = plotdir, height = 6, width = 8)

# cleanup --------
date()
closeAllConnections()
