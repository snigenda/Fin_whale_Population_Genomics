# Title: Read in the gathered results, generate summaries
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Feb 23 09:57:12 2021
# Modification: Use it for random starts dataset
# Modification: Allowing for partially missing data
# Date: Wed Aug 25 18:33:03 2021


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(ggplot2)

# def functions --------
split_bySFS <- function(statdt, nsim, nrep) {
    outlist = lapply(1:nsim, function(ii){
        temp = statdt[statdt$bootNum == ii,]
        if (nrow(temp) != nrep) {
            # stop('Missing correct replications per SFS. Please check.')
            print(paste0('bootSFS', ii, ':', nrow(temp), 'runs. Missing correct replications per SFS.'))
        }
        # sort the temp file
        temp = temp %>%
            dplyr::arrange(desc(MaxEstLhood))
        return(temp)
    })
    return(outlist)
}

# test convergence for a vector of parameter estimates (be more relaxed now, one fold of variation)
test_convergence <- function(paramval, cutoff = 1) {
    medval = median(paramval, na.rm = TRUE)
    rangeval = range(paramval, na.rm = TRUE)
    rangeval = rangeval[2]-rangeval[1]
    if (medval == 0 & rangeval == 0) {
        perval = 0
    } else {
        perval = rangeval/medval
    }
    if (perval < 0) {
        stop('Range should not be negative')
    }
    if (perval < cutoff & perval >= 0) {
        output = TRUE
    } else {
        output = FALSE
    }
    return(output)
}

plot_fscrep_res <- function(dt, varlist, prefix, bootid) {
    vars=c(varlist,"MaxEstLhood", "runNum")
    LL_data=dt[1,"MaxObsLhood"]
    newid=nrow(dt)
    forplot=dt[,vars] %>%
        dplyr::arrange(desc(MaxEstLhood)) %>%
        dplyr::mutate(new.id = 1:newid) %>%
        reshape2::melt(., id.vars="new.id")

    pp=ggplot(data=forplot, aes(x=new.id,y=value,color=variable)) +
        geom_point()+
        facet_wrap(.~ variable,ncol=1,scales="free_y")+
        geom_hline(data = subset(forplot, variable == "MaxEstLhood"), aes(yintercept = LL_data),color="red")+
        labs(title=paste(prefix, 'bootNum =', bootid))+
        theme_bw()
    return(pp)
}

plot_histogram <- function(bestlldt, estimate) {
    bestlldtp = bestlldt %>%
        reshape2::melt()
    estimatep = estimate %>%
        reshape2::melt()

    pp <- ggplot() +
        geom_histogram(data = bestlldtp, aes(x = value), color = 'black', fill = 'lightgrey') +
        facet_wrap(. ~ variable, scales = 'free') +
        geom_vline(data = estimatep, aes(xintercept = value), color = 'red') +
        theme_bw()
    return(pp)
}

get_quantile <- function(deltadt) {
    delta95 = data.frame(matrix(ncol = ncol(deltadt), nrow = 2))
    rownames(delta95) = c('low5','up95')
    colnames(delta95) = colnames(deltadt)
    for (ii in 1:ncol(deltadt)) {
        # d = quantile(deltadt[,ii], c(0.05,0.95), na.rm = TRUE)
        d = quantile(deltadt[,ii], c(0.025,0.975), na.rm = TRUE)
        delta95[,ii] = d
    }
    return(delta95)
}

get_CI <- function(delta95, estimate) {
    # check colnames match
    if (!identical(colnames(delta95), colnames(estimate))){
        stop('Wrong CI input colnames.')
    }
    est95 = data.frame(matrix(ncol = ncol(estimate), nrow = 2))
    rownames(est95) = c('low5','up95')
    colnames(est95) = colnames(estimate)
    for (ii in 1:ncol(delta95)) {
        est95['low5',ii] = estimate[1,ii] - delta95['up95',ii]
        est95['up95',ii] = estimate[1,ii] - delta95['low5',ii]
    }
    return(est95)
}
# def variables --------
# assume working directory from the called command
args = commandArgs(trailingOnly=TRUE)
infile = as.character(args[1]) # should be a file name not a path
# infile = '/Users/linmeixi/google_drive/finwhale/analyses/fsc26/param_btsp/neutral/resultsSummaries/1D.3Epoch.ENP/1D.3Epoch.ENP_v4_r10_bootFSC_n20rd_Summary_20210716.csv'
today = format(Sys.Date(), "%Y%m%d")
nsim = 100 # number of bootstrap SFS
nrep = 20 # number of replications per SFS

# models and parameters
models=c("1Epoch","2Epoch","3Epoch","3EpochTcur3","3EpochTcur2","4Epoch",
         "SplitNoMig","SplitSymMig","SplitAsyMig","SplitAsyMigTw2","SplitAsyMigTw3",
         "AncSplitAsyMig", "AncSplitIsoAsyMig", "AncSplitAsyMigChgGOC")

# note that here changed the order to be N + T + M for everything
# NPG = before the current GOC size
varlist_list = list(c("NE"),
                    c("NANC","NCUR","TCON"),
                    c("NCUR","NBOT","NANC","TBOT","TENDBOT"),
                    c("NCUR","NBOT","NANC","TENDBOT"),
                    c("NCUR","NBOT","NANC","TENDBOT"),
                    c("NCUR","NREC","NBOT","NANC","TBOT","TENDBOT","TANCBOT"),
                    c("NANC","NGOC","NENP","TDIV"),
                    c("NANC","NGOC","NENP","TDIV","MIG"),
                    c("NANC","NGOC","NENP","TDIV","N1M21","N2M12","MIG21","MIG12"),
                    c("NANC","NGOC","NENP","NBEC","TDIV","N1M21","N2M12","MIG21","MIG12"),
                    c("NANC","NGOC","NENP","NBEC","TDIV","N1M21","N2M12","MIG21","MIG12"),
                    c("NANC","NBOT","NGOC","NENP","TDIV","TA","N1M21","N2M12","MIG21","MIG12"),
                    c("NANC","NBOT","NGOC","NENP","TDIV","TA","TISO","N1M21","N2M12","MIG21","MIG12"),
                    c("NANC","NBOT","NGOC","NENP","NPG","TDIV","TA","TGULF","N1M21","N2M12","MIG21","MIG12"))

names(varlist_list)=models

sessionInfo()

# load data --------
# set the variables of interest
myprefix = stringr::str_split(infile, pattern = '_Summary')[[1]][1]
mymodel = stringr::str_split(myprefix, pattern = '\\.')[[1]][2]
# quit if not defined in varlist_list
if (! (mymodel %in% models)) {
    stop('Model not defined')
}
myvars = varlist_list[[mymodel]]
# load data
indt = read.csv(file = infile, stringsAsFactors = FALSE)
# split the dt for variants and output statistics
statdt = indt[-1,]
estimate = indt[1,]

# split the statdt into a list given nsim and nrep
statlist = split_bySFS(statdt, nsim, nrep)

# main --------
# look for convergence ========
# only look at the top 20 ranking runs (don't have more than 20 now)
converlist = lapply(statlist, function(dt) {
    outcon = unlist(lapply(myvars, function(var) {
        # print(var)
        output = test_convergence(dt[,var])
        if (output == FALSE) {
        bootid = dt[1, 'bootNum']
        print(paste0('WARNING: bootNum = ', bootid, ' bootSFS is not converging on ', var))
        }
        return(output)
    }))
    return(outcon)
})
converdt = do.call(rbind, converlist)
colnames(converdt) = myvars
write.csv(converdt,file = paste0(myprefix, '_SumConvergence_', today, '.csv'))

# plotting at random for three runs
bootids = sample(nsim,3)
pplist <- lapply(bootids, function(ii){
    pp = plot_fscrep_res(dt = statlist[[ii]], varlist = myvars, prefix = myprefix, bootid = ii)
    ggsave(pp, filename = paste0('plot_',myprefix, '_run', ii, '_',today,'.pdf'), height = 10, width = 8)
})

# pick the bestll data regardless of convergence ========
# NOTE statlist is already sorted from the bestLL
bestlldt = do.call(rbind, lapply(statlist, function(dt) {dt = dt[1,]}))
bestlldt = bestlldt %>%
    dplyr::mutate(DiffLhood = MaxObsLhood - MaxEstLhood)
write.csv(bestlldt,file = paste0(myprefix, '_SumBestLL_', today, '.csv'))

# get the 95% CI ========
# format the bestlldt and estimate dt
bestlldt = bestlldt[,c(myvars, 'DiffLhood')]
deltadt = data.frame(matrix(ncol = ncol(bestlldt), nrow = nrow(bestlldt)))
colnames(deltadt) = colnames(bestlldt)
estimate = estimate %>%
    dplyr::mutate(DiffLhood = MaxObsLhood - MaxEstLhood)
estimate = estimate[,c(myvars, 'DiffLhood')]

# plot a histogram for the bestlldt ########
pp <- plot_histogram(bestlldt, estimate)
ggsave(pp, filename = paste0('plotHIST_',myprefix, '_',today,'.pdf'), height = 9, width = 9)

# looping through each variable ########
for (ii in 1:ncol(deltadt)) {
    # check that colnames is the same
    if (colnames(bestlldt)[ii] == colnames(estimate)[ii]) {
        deltadt[,ii] = bestlldt[,ii] - estimate[1,ii]
    } else {
        stop('Wrong colnames for bestll and estimate')
    }
}

# get 95 CI
delta95 = get_quantile(deltadt)
est95 = get_CI(delta95, estimate)

# generate output dt ========
# OriginalM2SD = Original estimate minus 2 STDDEV
# OriginalP2SD = Original estimate plus 2 STDDEV
outmeasures = c('Model','Parameter','Original', 'Average', 'Min', 'Max',
    'QuantileLow025','QuantileHigh975',
    'delta95low25', 'delta95up975','SD',
    'OriginalM2SD', 'OriginalP2SD','low95CI', 'up95CI')
outdt = data.frame(matrix(nrow = ncol(bestlldt), ncol = length(outmeasures)))
colnames(outdt)=outmeasures

for (ii in 1:nrow(outdt)) {
    var = colnames(bestlldt)[ii]
    temp = bestlldt[,ii]
    temp = temp[!is.na(temp)]
    quanttemp = quantile(temp, c(0.025,0.975))
    outdt[ii, 'Model'] = myprefix
    outdt[ii, 'Parameter'] = var
    outdt[ii, 'Original'] = estimate[1,ii] # estimate should be a 1 row data frame
    outdt[ii, 'Average'] = mean(temp)
    outdt[ii, 'Min'] = min(temp)
    outdt[ii, 'Max'] = max(temp)
    outdt[ii, 'QuantileLow025'] = unname(quanttemp['2.5%'])
    outdt[ii, 'QuantileHigh975'] = unname(quanttemp['97.5%'])
    outdt[ii, 'delta95low25'] = delta95['low5',var]
    outdt[ii, 'delta95up975'] = delta95['up95',var]
    outdt[ii, 'SD'] = sd(temp)
    outdt[ii, 'OriginalM2SD'] = outdt[ii, 'Original'] - 2*outdt[ii, 'SD']
    outdt[ii, 'OriginalP2SD'] = outdt[ii, 'Original'] + 2*outdt[ii, 'SD']
    outdt[ii, 'low95CI'] = est95['low5',var]
    outdt[ii, 'up95CI'] = est95['up95',var]
}

# output files --------
print(outdt[,c('Model','Parameter','Original','QuantileLow025','QuantileHigh975','OriginalM2SD','OriginalP2SD','low95CI', 'up95CI')])
write.csv(outdt, file = paste0('fsc26_CI_', myprefix, '_',today,'.csv'))

# cleanup --------
date()
closeAllConnections()
