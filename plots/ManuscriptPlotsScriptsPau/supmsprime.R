#Libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape)
dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/Demography"
setwd(dir)

#Read msprime table runs summary 
msprime <- read.csv("msprime.csv", header = 1)
msprime <- melt(msprime)

#Plot
ggplot(msprime, aes(x=ï..Model, y= value, color= variable)) +
  geom_point(size=4) +
  facet_wrap(~factor(variable), scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("#A6E2DC","#97D49B","#E4A58F","#ECF87F")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.background =element_rect(fill="white"),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab(label = "") + ylab(label = "") 
