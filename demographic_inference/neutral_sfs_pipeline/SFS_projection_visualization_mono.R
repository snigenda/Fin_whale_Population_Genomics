##### Fin Whale Project
## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx
## NOTES: R 4.0,v1.0, July 30 2020
#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Join SFS projections for chromosomes and plot them
#--------------------------------------------------------------------------------------------------

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################


# Load the R packages
library(plyr)
#library(readr)
library(ggplot2)
library(RColorBrewer)
require(reshape2)
#library(ggpubr)

#Define functions ------------------------------ 
Monomorphic <- function(filenames, dataframe, pop){
  c <- 2
  for(file in filenames){
    mono <- read.table(file, header = T)
    #print(mono)
    dataframe <- cbind(dataframe,mono$HomREFcountPassingProjThreshold[mono$population == pop]) 
    colnames(dataframe)[c] <- paste0("VCF_",c)
    #print(dataframe)
    c <- c+1
  }
  return(dataframe)
}

Projection <- function(filenames, dataframe,fold){
  c <- 2
  for(file in filenames){
    proj <- read.table(file, header = T, skip = 1)
    proj2 <- melt(proj)
    #print(proj2[1:fold,2])
    dataframe <- cbind(dataframe,proj2[,2]) 
    colnames(dataframe)[c] <- paste0("VCF_",c)
    c <- c+1
  }
  return(dataframe)
}



#Define variables --------------------
dir <- "/u/project/rwayne/snigenda/finwhale/SFS/Neutral"
pops <- c("ENP","GOC")
dadidir <- "/u/project/rwayne/snigenda/finwhale/SFS/Neutral/SFS_projection/SNPs_easySFS_projection_01/dadi"

#Main---------------------

setwd(dir)

#Obtain Monomorphic sites files
mono.files <- list.files(paste0(dir,"/SFS_projection_Monomorphic"),pattern= ".perPop",full.names = T)
#print(mono.files)

#Obtain tables for each pop that contain projection values of each chromosome file
for(pop in pops){
  #Monomorphic
  mono <- data.frame("Frequency" = 0)
  mono <- Monomorphic(mono.files, mono,pop)
  #print(mono)
  mono$Counts <- rowSums(mono[c(2:97)])
  mono$Pop <- pop
  
  #EasySFS Results
  files <- scan(paste0(dir,"/SFS_projection_fsc/",pop,"_SFS_projection_files.txt"), what = character()) #Obtain list of complete paths of the files
  file1 <- read.table(files[1], skip = 1)  
  ssfold <- (length(file1)+1)/2 # Where is the fold? (remember there's the monomorphic column)
  projection.values <- data.frame("Frequency" = 0:(length(file1)-1)) #Obtain frequencies
  #print(projection.values)
  projection <- Projection(files,projection.values,ssfold) #Obtain all the values per file
  projection$Counts <- rowSums(projection[c(2:97)]) # Sum all the values per chromosome
  projection$Proportion <- projection$Counts/sum(projection$Counts) #For proportional SFS
  projection$Pop <- pop
  
  #Add monomorphic sites
  projection$Counts[1] <- projection$Counts[1] + mono$Counts
  
  assign(paste0("projection_",pop), projection)
  
  #Write new formats
  print(projection$Counts)
  counts <- projection$Counts
  sfs <- cbind.data.frame(split(counts, rep(1:length(file1), times=1)), stringsAsFactors=F)
  colnames(sfs) <- paste0("d0_",0:(length(file1)-1))
  ##Fastsimcoal
  sink(paste(dir,"/SFS_projection_fsc/",pop,"_MAFpop0.obs",sep=""))
  cat("1 observation\n")
  write.table(sfs,quote=F,row.names = F) # writes new sfs to the table
  sink()
  
  #dadi
  file <- list.files(dadidir, pattern = paste("^",pop,"-[0-9]+.sfs",sep=""), full.names = T)
  header <- readLines(file,n=1) # get first line of file
  sfs_dadi <- read.table(file,skip = 1,header = F)
  print(sfs_dadi)
  print(sfs_dadi[1,1:length(file1)])
  sfs_dadi[1,1:length(file1)] <- sfs 
  sink(paste(dir,"/SFS_projection_fsc/",pop,"-",as.character(length(sfs_dadi)-1),".sfs",sep=""))
  cat(header,"\n") # put the header back in
  write.table(sfs_dadi,quote=F,row.names = F,col.names = F) # writes new sfs to the table
  sink()
  totalSites <- sum(sfs_dadi[1,])
  sink(paste(dir,"/SFS_projection_fsc/",pop,"-",as.character(length(sfs_dadi)-1),".totalSiteCount.L.withMonomorphic.txt",sep=""))
  cat("pop\ttotalSites\n")
  cat(pop, totalSites,"\n")
  sink()
  
}


#Plot SFS
projection_allpops1 <- rbind(projection_ENP,projection_GOC)
projection_allpops2 <- projection_allpops1[-1,] 
  
colors <- brewer.pal(3,"Dark2")

pdf(paste0(dir,"/SFS_projection_fsc/","Final_fwhale_SFS.pdf"))

#Normal SFS
ggplot(projection_allpops2,aes(x=Frequency, y=Counts,fill=Pop)) + geom_bar(stat="identity")+ 
  theme_light() + scale_x_continuous(breaks=c(1:length(projection_ENP$Frequency-1)))+
  labs(title = "Folded SFS projection", x = "Alternative Allele count", y= "") + facet_grid(~Pop) +
  scale_fill_manual(values = colors) + guides(fill = F)

#Proportional SFS
ggplot(projection_allpops2,aes(x=Frequency, y=Proportion,fill=Pop)) + geom_bar(stat="identity")+ 
  theme_light() + scale_x_continuous(breaks=c(1:length(projection_ENP$Frequency-1)))+
  labs(title = "Folded proportional SFS ", x = "Alternative Allele count", y= "Proportion") + facet_grid(~Pop) +
  scale_fill_manual(values = colors) + guides(fill = F)

#Proportions 0 bin
ggplot(projection_allpops1,aes(x=Frequency, y=Counts,fill=Pop)) + geom_bar(stat="identity")+ 
  theme_light() + scale_x_continuous(breaks=c(0:length(projection_ENP$Frequency-1)))+
  labs(title = "Folded SFS projection", subtitle = "With monomorphic Sites",x = "Alternative Allele count", y= "") + facet_grid(~Pop) +
  scale_fill_manual(values = colors) + guides(fill = F)

dev.off()
