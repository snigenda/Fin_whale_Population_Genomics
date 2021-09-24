#Admixture

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(pophelper)
library(SNPRelate)
library(ggpubr)
library(gridExtra)


dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/PopStructure/ADMIXTURE/final"
gdsfile <- "../../JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds"
setwd(dir)
popmap = read.table("../../../popmap.txt", header = T)

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
colors <-c("#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F", "#BEBADA")
col <- c("#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")

afiles<- list.files(path=".", pattern = "*.Q", full.names=T)
alist <- readQ(files=afiles, indlabfromfile=TRUE)
tr1 <- tabulateQ(qlist=alist)

fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(alist,fn1))
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",sample.id)
lapply(alist, rownames)
pop <- popmap[c(3,2)]
pop <- pop

#Plot normal
p1 <- plotQ(alignK(sortQ(alist)),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11, showsp = F,
            clustercol=colors,legendtextsize=10,linesize=0.8,pointsize=3,
            showindlab=T,useindlab=T,
            grplab= pop, splab = spnames, grplabsize = 3,
            grplabheight = 5,subsetgrp=c("AK","BC","WA","OR","CA","GOC"))

grid.arrange(p1$plot[[1]])

#For each k
p1 <- plotQ(alignK(sortQ(alist[c(41:50)])),imgoutput="join",returnplot=T,
            exportplot=F,quiet=T, grplab= pop, grplabcol = "white", grplabsize = 1,
            clustercol = colors, showlegend = F, splabcol = "white",
            linesize=0.8,pointsize=3, subsetgrp=c("AK","BC","WA","OR","CA","GOC"),
            linecol = "white", pointcol = "white")
grid.arrange(p1$plot[[1]])


#Plot final
p1 <- plotQ(alignK(sortQ(alist[c(41:50)])),imgoutput="join",returnplot=T,
            exportplot=F,quiet=T, indlabsize = 8,  grplab= pop,
            grplabcol = "white", grplabsize = 1,
            clustercol = col, showlegend = F, splabcol = "white",
            linesize=0.8,pointsize=3, showindlab=T,useindlab=T, 
            subsetgrp=c("AK","BC","WA","OR","CA","GOC"),
            linecol = "white", pointcol = "white")
grid.arrange(p1$plot[[1]])
 #6*7

