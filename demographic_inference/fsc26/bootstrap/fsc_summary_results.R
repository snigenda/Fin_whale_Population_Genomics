# Title: Read in the summaried results, generate values to copy into word document
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Aug 29 11:22:31 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE,
        stringsAsFactors = FALSE)
setwd('~/Lab/finwhale_manuscript/')

library(dplyr)
library(stringr)

# def functions --------
# https://stackoverflow.com/questions/12688717/round-up-from-5
round2 = function(x, n=0) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z/10^n
    z*posneg
}

filter_parameters <- function(mydt) {
    model = str_split(mydt$Model[1], '\\.')[[1]][2]
    vars = varlist_list[[model]]
    varmatch = match(vars, mydt$Parameter)
    outdt = mydt[varmatch,]
    return(outdt)
}

convert_diploid <- function(mydt) {
    model = str_split(mydt$Model[1], '\\.')[[1]][2]
    vars = varlist_list[[model]] 
    vars = vars[str_detect(vars, '^N')]
    outdt = mydt
    for (ii in vars) {
        outdt[outdt$Parameter == ii, 3:5] = outdt[outdt$Parameter == ii, 3:5]/2
    }
    return(outdt)
}

get_interval <- function(mydt, time1, time2) {
    original = mydt[mydt$Parameter == time1, 'Original'] - mydt[mydt$Parameter == time2, 'Original']
    lowT = mydt[mydt$Parameter == time1, 'QuantileLow025'] - mydt[mydt$Parameter == time2, 'QuantileHigh975']
    highT = mydt[mydt$Parameter == time1, 'QuantileHigh975'] - mydt[mydt$Parameter == time2, 'QuantileLow025']
    return(c(original,lowT, highT))
}

time_3epoch <- function(mydt) {
    Tbot = get_interval(mydt, 'TENDBOT', 'TBOT')
    mydt[mydt$Parameter == 'TENDBOT',3:5] = Tbot
    return(mydt)
}

time_4epoch <- function(mydt) {
    Tbot = get_interval(mydt, 'TANCBOT', 'TENDBOT')
    Trec = get_interval(mydt, 'TENDBOT', 'TBOT')
    mydt[mydt$Parameter == 'TANCBOT',3:5] = Tbot
    mydt[mydt$Parameter == 'TENDBOT',3:5] = Trec
    return(mydt)
}

time_varfix <- function(mydt, var, tcur) {
    mydt[mydt$Parameter == var,3:5] = mydt[mydt$Parameter == var,3:5] - tcur
    return(mydt)
}

time_ancsplit <- function(mydt) {
    Ta = get_interval(mydt, 'TA', 'TDIV')
    mydt[mydt$Parameter == 'TA',3:5] = Ta
    return(mydt)
}

time_chggoc <- function(mydt) {
    Ta = get_interval(mydt, 'TA', 'TDIV')
    Td = get_interval(mydt,'TDIV', 'TGULF')
    mydt[mydt$Parameter == 'TA',3:5] = Ta
    mydt[mydt$Parameter == 'TDIV',3:5] = Td
    return(mydt)
}

time_anciso <- function(mydt) {
    Ta = get_interval(mydt, 'TA', 'TDIV')
    Td = get_interval(mydt,'TDIV', 'TISO')
    mydt[mydt$Parameter == 'TA',3:5] = Ta
    mydt[mydt$Parameter == 'TDIV',3:5] = Td
    return(mydt)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")

indir = './data/demography/fsc26/bootstrap/raw_data/'
outdir = './data/demography/fsc26/bootstrap/derived_data/'

modelorder <- as.character(c(1,7,2,8,5,9,20,21,6,10,12,13,14,16,15,17,18,19))
refiddt <- data.frame(refid = modelorder)

# parameters for each model
# models and parameters
# models and parameters
models=c("1Epoch","2Epoch","3Epoch","3EpochTcur3","3EpochTcur2","4Epoch",
         "SplitNoMig","SplitSymMig","SplitAsyMig","SplitAsyMigTw2","SplitAsyMigTw3",
         "AncSplitAsyMig", "AncSplitIsoAsyMig", "AncSplitAsyMigChgGOC")

# note that here changed the order to be N + T + M for everything
# NPG = before the current GOC size
varlist_list = list(c("NE"),
                    c("NANC","NCUR","TCON"),
                    c("NANC","NBOT","NCUR","TENDBOT","TBOT"),
                    c("NANC","NBOT","NCUR","TENDBOT"),
                    c("NANC","NBOT","NCUR","TENDBOT"),
                    c("NANC","NBOT","NREC","NCUR","TANCBOT","TENDBOT","TBOT"),
                    c("NANC","NENP","NGOC","TDIV"),
                    c("NANC","NENP","NGOC","TDIV","MIG"),
                    c("NANC","NENP","NGOC","TDIV","MIG21","MIG12"),
                    c("NANC","NBEC","NENP","NGOC","TDIV","MIG21","MIG12"), # NBEC = N_ENP NENP = N_ENP2
                    c("NANC","NBEC","NENP","NGOC","TDIV","MIG21","MIG12"),
                    c("NANC","NBOT","NENP","NGOC","TA","TDIV","MIG21","MIG12"), # NBOT=N_ANC2
                    c("NANC","NBOT","NENP","NGOC","TA","TDIV","TISO","MIG21","MIG12"),
                    c("NANC","NBOT","NENP","NPG","NGOC","TA","TDIV","TGULF","MIG21","MIG12")) # NPG = N_GOC NGOC =N_GOC2
names(varlist_list) = models

# load data --------
fsc_models <- read.csv(file = paste0('./scripts/demography/fsc26/bootstrap/fsc_output_list.csv')) %>%
    filter(!str_detect(refid,'#'))
fsc_models <- dplyr::left_join(refiddt, fsc_models, by = 'refid') %>%
    dplyr::mutate(modelname = paste(submodel, pop, sep = '.'))

# load the data frames ========
dtlist <- lapply(fsc_models$modelname, function(xx){
    myfile = list.files(path = paste0(indir, xx, '/'),
                        pattern = 'fsc26_CI_',
                        recursive = TRUE)
    print(myfile)
    mydt = read.csv(file = paste0(indir, xx, '/',myfile)) %>%
        dplyr::select(Model,Parameter,Original,QuantileLow025, QuantileHigh975)
    model = str_split(mydt$Model[1], '\\.')[[1]][2]
    mydt = filter_parameters(mydt)
    mydt = convert_diploid(mydt)
    if (model == "3Epoch") {mydt = time_3epoch(mydt)}
    if (model == "4Epoch") {mydt = time_4epoch(mydt)}
    if (model == "3EpochTcur3") {mydt = time_varfix(mydt, var = 'TENDBOT', tcur = 3)}
    if (model == "3EpochTcur2") {mydt = time_varfix(mydt, var = 'TENDBOT', tcur = 2)}
    if (model == "SplitAsyMigTw3") {mydt = time_varfix(mydt, var = 'TDIV', tcur = 3)}
    if (model == "SplitAsyMigTw2") {mydt = time_varfix(mydt, var = 'TDIV', tcur = 2)}
    if (model == "AncSplitAsyMig") {mydt = time_ancsplit(mydt)}
    if (model == "AncSplitIsoAsyMig") {mydt = time_anciso(mydt)}
    if (model == "AncSplitAsyMigChgGOC") {mydt = time_chggoc(mydt)}
    return(mydt)
})


# main --------

alldt = dplyr::bind_rows(dtlist)

# change print formats
alloutdt = alldt %>%
    dplyr::mutate(CI = ifelse(Parameter %in% c("MIG","MIG21","MIG12"),
                                  paste(formatC(QuantileLow025, format = "E", digits = 2), formatC(QuantileHigh975, format = "E", digits = 2), sep = " – "),
                                  paste(formatC(round2(QuantileLow025), format = "d"), formatC(round2(QuantileHigh975), format = "d"), sep = " – ")),
                  Estimate = ifelse(Parameter %in% c("MIG","MIG21","MIG12"),
                              formatC(Original, format = "E", digits = 2),
                              formatC(round2(Original), format = "d"))) %>%
    dplyr::select(Model, Parameter, Estimate, CI)

# output files --------
write.csv(x = alloutdt, file = paste0(outdir, 'MANUSCRIPTfsc_btsp_summary_n20rd_', today, '.csv'))

# cleanup --------
date()
closeAllConnections()
