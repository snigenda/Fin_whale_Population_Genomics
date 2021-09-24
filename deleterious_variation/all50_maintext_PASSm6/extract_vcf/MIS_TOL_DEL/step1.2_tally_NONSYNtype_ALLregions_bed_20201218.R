# Title: Tally up the NONSYN types
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Dec 14 17:07:52 2020

# preparation --------
rm(list = ls())
cat("\014")

options(echo = TRUE)

# def functions --------

# def variables --------
mywd = "/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/bedfiles"
setwd(mywd)

# main -------
dt = read.delim(file = 'JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_NONSYNtype.bed', header = FALSE)
# 'data.frame': 132525 obs. of  4 variables

colnames(dt) = c('CHROM', 'BEDSTART', 'END', 'NONSYNTYPE')
table(dt[,'NONSYNTYPE'], useNA = 'always')
# DELETERIOUS DELETERIOUS_(*WARNING!_Low_confidence) 
#       24922                                   8750 
#   TOLERATED                                   <NA> 
#       98750                                    103 

# output --------
dt = dt[!is.na(dt$NONSYNTYPE),]
dt_tol = dt[dt$NONSYNTYPE == 'TOLERATED',]
write.table(dt_tol, file = 'JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_TOL.bed',
	col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

dt_delall = dt[dt$NONSYNTYPE %in% c('DELETERIOUS','DELETERIOUS_(*WARNING!_Low_confidence)'),]
write.table(dt_delall, file = 'JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_DEL.bed',
	col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

# no warning
dt_delnow = dt[dt$NONSYNTYPE == 'DELETERIOUS',]
write.table(dt_delnow, file = 'JointCalls_all50_filterpassmiss_nonsyn_ALLregions_all_SIFT_DELnoW.bed',
	col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')