# Title: plot the projection preview output 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Jun  7 23:31:04 2020

# preparation --------
options(echo = FALSE)
rm(list = ls())
cat("\014")

require(ggplot2)
require(dplyr)

# def functions --------

# def variables --------
args = commandArgs(trailingOnly=TRUE)
data.dir = as.character(args[1])
plot.dir = as.character(args[2])
filePREFIX = as.character(args[3])

dir.create(plot.dir,recursive = T)
setwd(data.dir)

# load data --------
fileList=list.files(pattern=paste0(filePREFIX, ".easySFS.projPreview.R.format.txt"), path = data.dir, full.names = F)
# get populations (by stripping the )
pops = sub(paste0(".", filePREFIX, ".easySFS.projPreview.R.format.txt.*"), "", fileList)

# read in all data 
alldata = data.frame()
for (ii in seq_along(fileList)) {
    # read data 
    pop = pops[ii]
    tempdt = read.csv(file = fileList[ii]) %>% 
        dplyr::mutate(population = pop)
    alldata = rbind(alldata, tempdt)
}

# main --------
# get maxima 
maxima <- alldata %>%
  group_by(population) %>%
  filter(snps==max(snps))

print("The maxima for projections")
print(maxima)

# plot 
p0 <- ggplot(alldata,aes(x=projection,y=snps))+
  geom_point()+
  facet_wrap(~population,scales="free_x")+
  theme_bw()+
  scale_x_continuous(breaks=seq(2,56,2))
  
# want to find the maximum # of snps that still has good sample size 
ggsave(paste0(plot.dir, filePREFIX, ".easySFS.projections.allPops.pdf"),p0,height=5,width=10)

# cleanup --------
closeAllConnections()