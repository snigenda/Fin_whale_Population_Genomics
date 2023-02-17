#Admixture

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
# remotes::install_github('royfrancis/pophelper')
library(pophelper)
library(SNPRelate)
library(ggpubr)
library(gridExtra)


dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/PopStructure/ADMIXTURE/final"
gdsfile <- "../../JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds"
setwd(dir)
popmap = read.table("../../../popmap.txt", header = T)

# Modification: Save the source data
# Date: Sun Jan  8 16:43:42 2023
dir <- "/Users/linmeixi/Google Drive/My Drive/finwhale/analyses/PopStructure/all50/Minke/Admixture_20210318/maf10/"
gdsfile <- "~/Google Drive/My Drive/finwhale/analyses/PopStructure/all50/Minke/JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds"
popmap = read.csv(file = "/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/popmap_all50.csv", stringsAsFactors = F)
colnames(popmap) = c('sample','group','location')
setwd(dir)
###########


file <- snpgdsOpen(gdsfile)
#List of sample ids
sample.id <- read.gdsn(index.gdsn(file, "sample.id"))
#Geting pop and groups code
group_code <- popmap[,2]
pop_code <- popmap[,3]
#Close gds file
snpgdsClose(file)


colors = c("#D95F02","#1B9E77")

afiles<- list.files(path=".", pattern = "*.Q", full.names=T)
alist <- readQ(files=afiles, indlabfromfile=TRUE)
tr1 <- tabulateQ(qlist=alist)

fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(alist,fn1))
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",sample.id)
lapply(alist, rownames)
pop <- popmap[c(3,2)]


#Plot normal
p1 <- plotQ(alignK(alist[2]),returnplot=T,exportplot=F,
            # quiet=T,basesize=11, clustercol = colors,
            basesize=11, clustercol = colors,
            showyaxis =T, grplabsize=3,
            linesize=0.8,pointsize=3, showindlab=T,useindlab=T,
            grplab= pop, splab=spnames[2], legendtextsize=10,
            barbordersize = 0.001,
            barbordercolour="white",
            sortind="Cluster1",
            subsetgrp=c("AK","BC","WA","OR","CA","GOC"),
            showgrplab = T)
grid.arrange(p1$plot[[1]])



#Plot con subset ENP
p2 <- plotQMultiline(alignK(alist[2]),returnplot=T,exportplot=F,
            # quiet=T, barsize=0.98,grplab=pop[,2,drop=F],useindlab=T,
            barsize=0.98,grplab=pop[,2,drop=F],useindlab=T,
            sortind="Cluster1",clustercol = colors,
            subsetgrp=c("ENP"),
            basesize=30,grplabsize = 15, indlabsize = 12,
            barbordersize = 0.001,
            barbordercol="white",
            showyaxis = T)


grid.arrange(p2$plot[[1]][[1]])

p1 <- plotQ(alignK(alist[2]),returnplot=T,exportplot=F,
            # quiet=T,basesize=0, clustercol = colors,
            basesize=0, clustercol = colors,
            showyaxis =T, grplabsize=3,
            linesize=0.8,pointsize=3, showindlab=T,useindlab=T,
            grplab= pop[,1,drop=F], splab="",
            barbordersize = 0.001,
            barbordercolour="white",
            sortind="Cluster1",
            subsetgrp=c("AK","BC","WA","OR","CA"),
            showgrplab = T, indlabsize = 12)
grid.arrange(p1$plot[[1]])


#Plot con subset GOC
p2 <- plotQMultiline(alignK(alist[2]),returnplot=T,exportplot=F,
                     # quiet=T, barsize=0.98,grplab=pop[,2,drop=F],useindlab=T,
                     barsize=0.98,grplab=pop[,2,drop=F],useindlab=T,
                     sortind="Cluster1",clustercol = colors,
                     subsetgrp=c("GOC"),
                     basesize=30,grplabsize = 15, indlabsize = 12,
                     barbordersize = 0.001,
                     barbordercol="white",
                     grplabbgcol = "white")



grid.arrange(p2$plot[[1]][[1]])


##4*6


# Modification: Generate source data
# Date: Mon Jan  9 14:18:50 2023
# In the final data, we used one example run the `all50_pass_maf10.K2.iter10.Q`
fig1cdt = alignK(alist[2])[[1]]
write.csv(fig1cdt, file = '~/Lab/fin_whale/FinWhale_PopGenomics_2021/source_data/Fig1c.csv')





