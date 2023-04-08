##### Moreno Lab
## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx
## NOTES: R 4.0,v1.0, July 23 2020
#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Put together the preview of SFS projections from each chromosomes/scaffold
#--------------------------------------------------------------------------------------------------

# Usage in hoffman: Rscript SFS_preview_v1.R

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################


# Load the R packages
library(plyr)
library(readr)
library(ggplot2)
library(RColorBrewer)

#Define functions ------------------------------ 
readtables <- function(filenames, dataframe){
  c <- 1
  for(file in filenames){
    dataframe <- cbind(dataframe,read_csv(file,skip = 1,col_names = c("prj", paste0("SNPs_", c)) )[2]) 
    c <- c+1
  }
  return(dataframe)
}


#Define variables --------------------
dir <- "/u/project/rwayne/snigenda/finwhale/SFS/Neutral/SNPs"
dir.create("/u/project/rwayne/snigenda/finwhale/SFS/Neutral/Preview")
outputdir <- "/u/project/rwayne/snigenda/finwhale/SFS/Neutral/Preview"
pops <- c("ENP","GOC")


#Main---------------------

#Obtain tables for each pop that contain projection values of each chromosome file
for(pop in pops){
  preview.files <- list.files(dir,pattern= pop,full.names = TRUE)
  preview.values <- read_csv(preview.files[1], col_names = c("Projection","snp"), skip = 1)[1]
  preview <- readtables(preview.files,preview.values)
  write.table(preview, paste0(outputdir,"/SFS_preview_",pop,"_Sergio.txt"),quote = F, row.names = F)
  preview$SNPs <- rowSums(preview[c(2:97)])
  preview$pop <- pop
  assign(paste0("preview_",pop), preview)
}

#Join populations
preview_allpops <- rbind(preview_ENP, preview_GOC)
colors <- brewer.pal(3,"Dark2")

#Plot preview
p1 <- ggplot(preview_allpops,aes(x = Projection, y=SNPs), group = 1, color=pop) + geom_line() + geom_point() + 
  facet_wrap(~pop,scales="free_x") + theme_light() +  scale_x_continuous(breaks=seq(2,80,2)) +
  labs(title = "Site Frequency Spectrum Projection Preview", y="Number of segregating sites") + 
  guides(color = F) + scale_color_manual(values = colors)

ggsave(paste0(outputdir,"/FinWhale_easySFS_projections_allPops.pdf"),p1,height=5,width=10)

