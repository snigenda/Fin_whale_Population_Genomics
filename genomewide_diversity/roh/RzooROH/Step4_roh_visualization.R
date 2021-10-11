##### Moreno Lab
## Author: Paulina G. Nuñez Valencia, LCG-UNAM, Mexico: pnunez@lcg.unam.mx
## NOTES: R 4.0,v1.0, July 4th, 2020
#--------------------------------------------------------------------------------------------------
#                                            GOAL
# Plot rsults for ROH of mix10R model
#--------------------------------------------------------------------------------------------------

#################################################################################################
#                                        MAIN PROGRAM                                           #
#################################################################################################

#Loading packges 
require(RZooRoH)
require(R.utils)
require(ggplot2)
require(RColorBrewer)
require(reshape)
library(viridis)



#Set working directory
setwd("/u/home/p/pnunez/project-rwayne/ROHs/RZOOROH")
source("zooplot_prophbd_v1.R")
source("zooplot_reverse_v2.R")
source("zooplot_individuals_v1.R")

#Colors
pops <- c("ENP","GOC") 
model <- "mix10R_base3"

### load in each population's results and assign to "pop_results":
for(pop in pops){
  load(paste("./",model,"/",pop,"_ZooROH_Results.RData",sep="")) 
  assign(paste(pop,"_results",sep=""), results)
}

popResultsList=list(ENP=ENP_results,GOC=GOC_results)

pdf(paste0(model,"_pops.pdf"))
zooplot_prophbd_v1(popResultsList,style='barplot')
#png(paste0(model,"_pops_lines.png"))
#zooplot_prophbd_v1(popResultsList,style='lines', cumulative = TRUE)
#dev.off()
#png(paste0(model,"_pops_reverse.png"))
zooplot_partitioning_v1(popResultsList, plotids = FALSE,ylim=c(0,0.5), nonhbd = FALSE)
#dev.off()
#png(paste0(model,"_pops_chr1.png"))
zooplot_hbdseg(popResultsList, chr=1,coord=c(10000000,50000000))
dev.off()

colors <- brewer.pal(n = 4, "Dark2")
grop <- function(group, color){
  pdf(paste0(model,"_",group,".pdf"))
  allseg <- rohbd(zres = popResultsList[[group]])
  #write.table(allseg,file=paste0(model,"_",group,"segments.txt"),quote=F, row.names=F, sep="\t")
  p1 <- ggplot(allseg, aes(x= length)) + geom_histogram(fill = color) + theme_light() +
    labs(title = paste0("Length of HBD segments in ", group, " group"), x = "Length of HBD segment ", y = "Frequency")
  print(p1)
  #dev.off()
  #png(paste0(model,"_",group,"_boxplot.png"))
  zooplot_prophbd_v1(popResultsList[[group]], cols = color, style = 'boxplot')
  #dev.off()
  #png(paste0(model,"_",group,"_cumulative.png"))
  zooplot_individuals_v1(popResultsList[[group]], cumulative = TRUE, col = color)
  #dev.off()
  #png(paste0(model,"_",group,"_reverse.png"))
  print(zooplot_partitioning_v1(popResultsList[[group]], ylim = c(0,0.5), nonhbd = FALSE))
  dev.off()
}

grop("ENP", colors[1])
grop("GOC", colors[2])

