library(ggplot2)
library(dplyr)
library(ggpubr)
dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/Demography"

setwd(dir)

#################

load_dadi_sfs <- function(dadisfs, pop, id) {
  sfs <- read.table(dadisfs,skip = 1,header = FALSE, colClasses = 'numeric')
  rownames(sfs) = c('count', 'mask')
  nsites = dim(sfs)[2]-1
  sfs = t(sfs) %>%
    as.data.frame() %>%
    dplyr::mutate(site_freq = seq(0,nsites),
                  sfsid = id, sfspop=pop) %>%
    dplyr::filter(mask == 0)
  sfssum = sum(sfs$count)
  sfs = sfs %>%
    dplyr::mutate(prop = count/sfssum)
  return(sfs)
}

plot_fit_sfs_color <- function(sfsplot, prefix) {
  pp = ggplot(data = sfsplot, aes(x = site_freq, y = prop)) +
    geom_bar(stat="identity", position=position_dodge(),  aes(alpha = factor(sfsid), fill= sfspop))+ 
    scale_alpha_manual(values = c( 1,0.7,0.5))+
    theme_bw()+
    xlab("") +
    ylab("")  +
    scale_fill_manual(values = c("#1B9E77")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.background =element_rect(fill="white"),
          legend.position = "none")
  return(pp)
}

##ENP
ENPSFSdata <- load_dadi_sfs("ENP-44.sfs","ENP"," Data")
ENPSFSmodel <- load_dadi_sfs("ENP.best.sfs","ENP","dadi best model")
ENPSFSmsprime <- load_dadi_sfs("ENP.msprime.sfs","ENP","msprime simulation")
ENP <- rbind(ENPSFSdata,ENPSFSmodel,ENPSFSmsprime)
plot_fit_sfs_color(ENP, "")


##8*4
pp = ggplot(data = ENP, aes(x = site_freq, y = prop)) +
  geom_bar(stat="identity", position=position_dodge(),  aes(alpha = factor(sfsid), fill= sfspop))+ 
  scale_alpha_manual(values = c( 1,0.7,0.5))+
  theme_bw()+
  xlab("") +
  ylab("")  +
  scale_fill_manual(values = c("#1B9E77")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.background =element_rect(fill="white")) +
  guides(fill="none") + 
  labs(alpha="")

as_ggplot(get_legend(pp))

##4*4
