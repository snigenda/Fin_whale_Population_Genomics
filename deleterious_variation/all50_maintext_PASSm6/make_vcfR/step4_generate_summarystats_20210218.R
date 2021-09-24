# Title: Generate summary statistics for both mutdt and typedt for the newly filtered dataset
# NOTE: All of the typedt and mutdt here are 44 samples. Only variants that are segregating in these 44 samples were counted
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Feb 19 11:57:24 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(reshape2)
library(ggplot2)

# def functions --------
# get summary to report
get_summary_typedt <- function(dt) {
    # for everything, get a summary
    longdt = reshape2::melt(dt, id.vars = c("SampleId", "PopId", "SubPopId", "EffectType"))
    longdt = longdt %>%
        dplyr::group_by(PopId,EffectType, variable) %>%
        dplyr::summarise(meanVal = mean(value),
                         medianVal = median(value))
    # cast by popId
    castdt1 = longdt %>%
        reshape2::dcast(formula =  EffectType + variable ~ PopId, value.var = 'meanVal') %>%
        dplyr::mutate(ENPmGOC = ENP - GOC) %>%
        dplyr::mutate(ENPmGOCdENP = ENPmGOC/ENP,
                      ENPmGOCdGOC = ENPmGOC/GOC)
    castdt2 = longdt %>%
        reshape2::dcast(formula =  EffectType + variable ~ PopId, value.var = 'medianVal') %>%
        dplyr::mutate(ENPmGOC = ENP - GOC) %>%
        dplyr::mutate(ENPmGOCdENP = ENPmGOC/ENP,
                      ENPmGOCdGOC = ENPmGOC/GOC)
    return(list(castdt1, castdt2))
}

get_summary_nopop_typedt <- function(dt) {
    # for everything, get a summary
    longdt = reshape2::melt(dt, id.vars = c("SampleId", "PopId", "SubPopId", "EffectType"))
    longdt = longdt %>%
        dplyr::group_by(EffectType, variable) %>%
        dplyr::summarise(meanVal = mean(value),
                         medianVal = median(value))
    return(longdt)
}

# get summary to report for mutation types
get_summary_mutdt <- function(dt) {
    # for everything, get a summary
    longdt = reshape2::melt(dt, id.vars = c("SampleId", "PopId", "SubPopId", "MutType"))
    longdt = longdt %>%
        dplyr::group_by(PopId,MutType, variable) %>%
        dplyr::summarise(meanVal = mean(value),
                         medianVal = median(value))
    # cast by popId
    castdt1 = longdt %>%
        reshape2::dcast(formula =  MutType + variable ~ PopId, value.var = 'meanVal') %>%
        dplyr::mutate(ENPmGOC = ENP - GOC) %>%
        dplyr::mutate(ENPmGOCdENP = ENPmGOC/ENP,
                      ENPmGOCdGOC = ENPmGOC/GOC)
    castdt2 = longdt %>%
        reshape2::dcast(formula =  MutType + variable ~ PopId, value.var = 'medianVal') %>%
        dplyr::mutate(ENPmGOC = ENP - GOC) %>%
        dplyr::mutate(ENPmGOCdENP = ENPmGOC/ENP,
                      ENPmGOCdGOC = ENPmGOC/GOC)
    return(list(castdt1, castdt2))
}

get_summary_nopop_mutdt <- function(dt) {
    # for everything, get a summary
    longdt = reshape2::melt(dt, id.vars = c("SampleId", "PopId", "SubPopId", "MutType"))
    longdt = longdt %>%
        dplyr::group_by(MutType, variable) %>%
        dplyr::summarise(meanVal = mean(value),
                         medianVal = median(value))
    return(longdt)
}

# get mean normalized variations (by mean value)
get_mean_rawdt <- function(rawdt) {
    # from nopopdt, get Average CalledCount and CalledAlleleCount
    meandt = rawdt %>%
        dplyr::select(CalledCount, MutType) %>%
        group_by(MutType) %>%
        summarise(meanCalledCount = mean(CalledCount),
                  meanCalledAllele = 2*meanCalledCount)
    outdt = dplyr::left_join(x=rawdt, y = meandt, by = 'MutType') %>%
        dplyr::mutate(normHomRef = HomRefPer * meanCalledCount,
                      normHomAlt = HomAltPer * meanCalledCount,
                      normHet = HetPer * meanCalledCount,
                      normRefAllele = RefAllelePer * meanCalledAllele,
                      normAltAllele = AltAllelePer * meanCalledAllele)
    return(outdt)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
cdstype = 'ALLregions'
mutprefix = c("syn", "nonsyn", "nonsynTOL", "nonsynDEL", "LOF")
typeprefix = c("BenignRM", "Damaging")

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR/', dataset, ref, sep = '/')
setwd(workdir)
outdir = './derive_data/statssum_table/'
dir.create(outdir)

sessionInfo()

# load data --------
rawdt = readRDS(file = './derive_data/sum_table/all50_Minke_ALLregions_SUMtable_Seg_PASSm6_20210218.rds')

# main --------
# convert percentages to number of average variations ========
normdt = get_mean_rawdt(rawdt)
# output this
saveRDS(normdt, file = './derive_data/sum_table/all50_Minke_ALLregions_SUMtable_Seg_PASSm6_Norm_20210219.rds')
# subset to two groups
mutdt = normdt %>%
    dplyr::filter(MutType %in% mutprefix)
typedt = normdt %>%
    dplyr::filter(MutType %in% typeprefix) %>%
    dplyr::rename(EffectType = MutType)

# get and output summary ========
typedtsum = get_summary_typedt(typedt)
write.csv(typedtsum[[1]], file = paste0(outdir, "meandt_EffectType_", today, ".csv"))
write.csv(typedtsum[[2]], file = paste0(outdir, "mediandt_EffectType_", today, ".csv"))

mutdtsum = get_summary_mutdt(mutdt)
write.csv(mutdtsum[[1]], file = paste0(outdir, "meandt_MutType_", today, ".csv"))
write.csv(mutdtsum[[2]], file = paste0(outdir, "mediandt_MutType_", today, ".csv"))

# get and output summary for all populations ========
typenopop = get_summary_nopop_typedt(typedt)
mutnopop = get_summary_nopop_mutdt(mutdt)
write.csv(typenopop, file = paste0(outdir, "nopopdt_EffectType_", today, ".csv"))
write.csv(mutnopop, file = paste0(outdir, "nopopdt_MutType_", today, ".csv"))

# some important values ========
imp_typedtsum = lapply(typedtsum, function(xx) {
    yy = xx[stringr::str_detect(xx$variable, 'Per$')|stringr::str_detect(xx$variable, '^norm'),]
})
write.csv(imp_typedtsum[[1]], file = paste0(outdir, "IMPmeandt_EffectType_", today, ".csv"))
write.csv(imp_typedtsum[[2]], file = paste0(outdir, "IMPmediandt_EffectType_", today, ".csv"))

imp_mutdtsum = lapply(mutdtsum, function(xx) {
    yy = xx[stringr::str_detect(xx$variable, 'Per$')|stringr::str_detect(xx$variable, '^norm'),]
})
write.csv(imp_mutdtsum[[1]], file = paste0(outdir, "IMPmeandt_MutType_", today, ".csv"))
write.csv(imp_mutdtsum[[2]], file = paste0(outdir, "IMPmediandt_MutType_", today, ".csv"))

# cleanup --------
closeAllConnections()


