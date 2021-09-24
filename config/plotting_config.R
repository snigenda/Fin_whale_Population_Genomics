# Title: configuration for plotting
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jun 12 11:55:50 2020

# def functions --------

# def variables --------
# plotting colors
mypops = c("ENP", "GOC")
# mycolors = RColorBrewer::brewer.pal(n = 3, name = "Dark2")[1:2]
mycolors = c("#1B9E77","#D95F02")
names(mycolors) = mypops

mylocs = c("AK", "BC", "CA", "OR", "WA", "GOC")
# loccolors = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[c(3:7,2)]
loccolors = c("#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#D95F02")
names(loccolors) = mylocs
subpoporder = c("AK", "BC", "WA", "OR", "CA", "GOC") # order for subpopulations 

myspecies = c("BalAcu", "BalMus", "EubGla", "MegNov", "BalPhy")
# speccolors = pals::alphabet(n=26)[c(1,2,4,10,9)]
speccolors = c("#F0A0FF","#0075DC","#4C005C","#94FFB5","#808080")
names(speccolors) = myspecies

# pals::pal.bands(speccolors)

# discarding individuals
# ENPOR12: low genotype depth
# ENPCA01: Admixture individuals
# ENPCA09: Admixture individuals
# GOC010: Admixture individuals
# GOC080: Related to GOC006 in relateness analyses
# GOC111: Related to GOC025 in relateness analyses
discardind = c("ENPOR12", "ENPCA01", "ENPCA09", "GOC010", "GOC080", "GOC111")
