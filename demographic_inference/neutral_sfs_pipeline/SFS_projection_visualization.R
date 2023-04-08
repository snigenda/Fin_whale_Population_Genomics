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
library(ggplot2)
library(RColorBrewer)
require(reshape2)
#library(ggpubr)

#Define functions ------------------------------ 
readtables <- function(filenames, dataframe,fold){
  c <- 1
  for(file in filenames){
    proj <- read.table(file, header = T, skip = 1)
    proj2 <- melt(proj)
    dataframe <- cbind(dataframe,proj2[2:fold,2]) 
    colnames(dataframe)[c+1] <- paste0("VCF_",c)
    c <- c+1
  }
  return(dataframe)
}

#Define variables --------------------
dir <- "/u/project/rwayne/snigenda/finwhale/SFS/Neutral/SFS_projection_fsc"
pops <- c("ENP","GOC")


#Main---------------------

setwd(dir)

#Obtain tables for each pop that contain projection values of each chromosome file
for(pop in pops){
  files <- scan(paste0(pop,"_SFS_projection_files.txt"), what = character()) #Obtain list of complete paths of the files
  file1 <- read.table(files[1], skip = 1)  
  ssfold <- (length(file1)+1)/2 # Where is the fold? (remember there's the monomorphic column)
  projection.values <- data.frame("Frequency" = 2:ssfold-1) #Obtain frequencies
  projection <- readtables(files,projection.values,ssfold) #Obtain all the values per file
  projection$Counts <- rowSums(projection[c(2:97)]) # Sum all the values per chromosome
  write.table(projection, paste0("SFS_projection_",pop),quote = F, row.names = F) #Check point
  projection$Proportion <- projection$Counts/sum(projection$Counts) #For proportional SFS
  projection$Pop <- pop
  assign(paste0("projection_",pop), projection)
}

#Join populations
projection_allpops <- rbind(projection_ENP, projection_GOC)
colors <- brewer.pal(3,"Dark2")

pdf("SFS.pdf")

#Normal SFS
ggplot(projection_allpops,aes(x=Frequency, y=Counts,fill=Pop)) + geom_bar(stat="identity")+ 
  theme_light() + scale_x_continuous(breaks=c(1:length(projection_ENP$Frequency)))+
  labs(title = "Folded SFS projection", x = "Alternative Allele count", y= "") + facet_grid(~Pop) +
  scale_fill_manual(values = colors) + guides(fill = F)

#Proportional SFS
ggplot(projection_allpops,aes(x=Frequency, y=Proportion,fill=Pop)) + geom_bar(stat="identity")+ 
  theme_light() + scale_x_continuous(breaks=c(1:length(projection_ENP$Frequency)))+
  labs(title = "Folded proportional SFS ", x = "Alternative Allele count", y= "Proportion") + facet_grid(~Pop) +
  scale_fill_manual(values = colors) + guides(fill = F)

dev.off()
