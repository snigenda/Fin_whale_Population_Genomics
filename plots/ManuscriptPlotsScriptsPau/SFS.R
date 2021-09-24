#Libraries
library(ggplot2)
library(dplyr)

dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/Demography"

setwd(dir)

#################

#Load sfs dadi
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

#Plot function
plot_fit_sfs_color <- function(sfsplot, prefix) {
  pp = ggplot(data = sfsplot, aes(x = site_freq, y = prop)) +
    geom_bar(stat="identity", position=position_dodge(),  aes(alpha = factor(sfsid), fill= sfspop))+ 
    scale_alpha_manual(values = c( 1,0.3))+
    facet_grid(. ~ sfspop, scales = "free_x")+
    theme_bw()+
    xlab("") +
    ylab("")  +
    scale_fill_manual(values = c("#1B9E77","#D95F02")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.background =element_rect(fill="white"),
          legend.position = "none")
  return(pp)
}


#Plot by epoch
plot_fit_sfs_color_epoch <- function(sfsplot, prefix) {
  pp = ggplot(data = sfsplot, aes(x = site_freq, y = prop, fill=sfsid)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw()+
    xlab("") +
    ylab("")  +
    scale_fill_manual(values = c("#1B9E77","#5ab4ac","#c7eae5","#d8b365","#f6e8c3")) +
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
ENPSFSmodel <- load_dadi_sfs("ENP.best.sfs","ENP","3 Epoch")

###Best fit
ENP <- rbind(ENPSFSdata,ENPSFSmodel)
plot_fit_sfs_color(ENP, "")

##All epochs
ENPSFSmodel1E <- load_dadi_sfs("ENP.1Epoch.sfs","ENP","1 Epoch")
ENPSFSmodel2E <- load_dadi_sfs("ENP.2Epoch.sfs","ENP","2 Epoch")
ENPSFSmodel4E <- load_dadi_sfs("ENP.4Epoch.sfs","ENP","4 Epoch")

ENP <- rbind(ENPSFSdata,ENPSFSmodel,ENPSFSmodel1E,ENPSFSmodel2E,ENPSFSmodel4E)

plot_fit_sfs_color_epoch(ENP,"")


##GOC
GOCSFSdata <- load_dadi_sfs("GOC-30.sfs","GOC","Data")
GOCSFSmodel <- load_dadi_sfs("GOC.best.sfs","GOC","Model")

###Best fit
GOC <- rbind(GOCSFSdata,GOCSFSmodel)
plot_fit_sfs_color(GOC, "")


##Both best fit
sfsplot <- rbind(ENPSFSdata,ENPSFSmodel,GOCSFSdata,GOCSFSmodel)
plot_fit_sfs_color(sfsplot, "")


##4.5*6

####Legend
pp = ggplot(data = sfsplot, aes(x = site_freq, y = prop)) +
  geom_bar(stat="identity", position=position_dodge(),  aes(alpha = factor(sfsid), fill= sfspop))+ 
  scale_alpha_manual(values = c( 1,0.3))+
  facet_grid(. ~ sfspop,scales="free_x")+
  theme_bw()+
  xlab("Derived Allele Count") +
  ylab("Proportion")  +
  scale_fill_manual(values = c("#1B9E77","#D95F02")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        strip.background =element_rect(fill="white")) +
  guides(fill=F) +
  labs(alpha= "")

as_ggplot(get_legend(pp))

##

pp = ggplot(data = ENP, aes(x = site_freq, y = prop, fill=sfsid)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw()+
  xlab("") +
  ylab("")  +
  scale_fill_manual(values = c("#1B9E77","#5ab4ac","#c7eae5","#d8b365","#f6e8c3")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.background =element_rect(fill="white")) +
  labs(fill="")

as_ggplot(get_legend(pp))

##4*4