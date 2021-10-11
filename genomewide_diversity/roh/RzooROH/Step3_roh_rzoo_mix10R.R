##### Runs of homozigosity 
## NOTES: R,v3.6, June 22, 2020
## Author: Jacqueline Robinson, modified for fin whale data by Paulina Nunez Valencia (pnunez@lcg.unam.mx) 
#==================================================================================================
#                                            GOAL
# Obtain runs of homozigosity for each individual for each pop
#--------------------------------------------------------------------------------------------------
# Parameters:
# 1) -dir complete path where the files to work with are
# Ouput:
# 1) R.data with Rzoo model for each pop
# 2) pdf with plots
# Execution:
# Rscript Rzoo_mix10R.R /u/home/p/pnunez/project-rwayne/ROHs/RZOOROH
#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################

# Load the R packages
library(RZooRoH)
library(R.utils)

#Define variables --------------------
args <- commandArgs(trailingOnly=TRUE)
dir <- as.character(args[1])

#Main  --------------------------
setwd(dir)
modelName <- "mix10R"
pops=c("GOC") 

for(pop in pops){
  print(paste("Starting pop",pop))
  out.dir <- paste0(dir,"/",modelName,"/")
  dir.create(out.dir,showWarnings = F)
  
  zooInput <- zoodata(paste0(dir,"/",pop,"/JointCalls_08_B_VariantPASS_",pop,"_allcontig.gen.gz"), zformat="gp",samplefile = paste0(dir,"/",pop,"/",pop,"_samples.txt"),min_maf =0.05)  
  model <- zoomodel() 
  model

  results <- zoorun(model, zooInput) 
  save(results, file = paste0(out.dir,pop,"_ZooROH_Results.RData"))

}
