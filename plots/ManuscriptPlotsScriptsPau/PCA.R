#PCA

library(dplyr)
library(ggplot2)
library(SNPRelate)
library(ggpubr)
library(ggfortify)
library(ggrepel)


dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/PopStructure"

gds <- "./JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds"

# Modification: Save the source data
# Date: Sun Jan  8 16:43:42 2023
dir <- "~/Lab/fin_whale/FinWhale_PopGenomics_2021/source_data/"
gds <- "~/Google Drive/My Drive/finwhale/analyses/PopStructure/all50/Minke/JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10.gds"

setwd(dir)

format_pca <- function(pca, popmap, poplevel) {
  pc.percent <- pca$varprop*100
  tab <- data.frame(sample.id = pca$sample.id,
                    loc = popmap[match(pca$sample.id, popmap$sample), poplevel],
                    pop = popmap[match(pca$sample.id, popmap$sample), "group"],
                    PC1 = pca$eigenvect[,1], # the first eigenvector
                    PC2 = pca$eigenvect[,2], # the second eigenvector
                    PC3 = pca$eigenvect[,3], # the third eigenvector
                    stringsAsFactors = FALSE)
  lbls <- paste0("PC", 1:4, "(", format(pc.percent[1:4], digits=2), "%", ")")
  return(list(tab, lbls))
}

genofile <- snpgdsOpen(gds, readonly = TRUE)
# popmap = read.table("../popmap.txt", header = T)
popmap = read.csv(file = "<homedir>/fin_whale/scripts_analyses/config/popmap_all50.csv", stringsAsFactors = F)
colnames(popmap) = c('sample','group','location')

pcaout = snpgdsPCA(genofile, num.thread=8, autosome.only=FALSE, verbose = FALSE)
pcaformatted = format_pca(pcaout, popmap = popmap, poplevel = "location")
tab = pcaformatted[[1]];lbls = pcaformatted[[2]]

loccolors = c("#7570B3","#E7298A","#66A61E","#D95F02","#E6AB02","#A6761D")
mycolors <-c("#1B9E77","#D95F02")
tab$Fac <- "PCA, all locations"


ggplot(tab, aes(x=PC1, y =PC2, color= loc)) +
  geom_point(size = 5) +
  labs(x = lbls[1], y = lbls[2]) +
  scale_color_manual(values = loccolors) +
  theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())


tab2 <- tab[c(24,33),]


ggplot(tab, aes(x=PC1, y =PC2, color= loc)) +
  geom_point(size = 5) +
  labs(x = lbls[1], y = lbls[2]) +
  scale_color_manual(values = loccolors) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   legend.position = "none") +
  geom_label_repel(data=tab2, aes(label = sample.id),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')


#Legend

loccolors2 = c("#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#D95F02")

p2 <- ggplot(tab, aes(x=PC1, y =PC2, color= loc)) +
  geom_point(size = 5) +
  labs(x = lbls[1], y = lbls[2]) +
  scale_color_manual(values = loccolors2, breaks = c("AK","BC","CA","OR","WA","GOC")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(color="")

as_ggplot(get_legend(p2))

## 5*5

closefn.gds(genofile)

# Modification: Meixi Lin
# Date: Sun Jan  8 16:41:09 2023
# Save the source data
write.csv(tab, file = 'Fig1b.csv')







