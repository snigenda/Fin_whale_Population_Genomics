library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/Demography"

setwd(dir)

#################

load_dadi_sfs <- function(dadisfs, pop, id) {
  sfs <- read.table(dadisfs,skip = 1,header = FALSE, colClasses = 'numeric')
  rownames(sfs) = c('value', 'variable')
  nsites = dim(sfs)[2]-1
  sfs = t(sfs) %>%
    as.data.frame() %>%
    dplyr::mutate(site_freq = seq(0,nsites),
                  sfsid = id, sfspop=pop, soft="dadi") %>%
    dplyr::filter(variable == 0)
  sfssum = sum(sfs$value)
  sfs = sfs %>%
    dplyr::mutate(prop = value/sfssum)
  return(sfs)
}


load_fsc_sfs <- function(fscsfs,pop, id,read){
  if(read == 1){
    sfs <- read.table(fscsfs, header = T, skip = 1)
  }
  else{
    sfs <- read.table(fscsfs, header = T)
  }
  sfs <- melt(sfs)
  sfs <- sfs[c(2:((length(sfs[,1])+1)/2)),]
  sfs$prop <- sfs$value/sum(sfs$value) #For proportional SFS
  sfs$sfsid <- id
  sfs$sfspop <- pop
  sfs$site_freq <- as.numeric(sapply(strsplit(as.character(sfs$variable),"_"),"[",2))
  sfs$soft <- "fastsimcoal2"
  return(sfs)
}


plot_fit_sfs_color_epoch <- function(sfsplot, prefix, colors) {
  pp = ggplot(data = sfsplot, aes(x = site_freq, y = prop, fill=sfsid)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw()+ facet_grid(vars(sfspop), vars(soft)) +
    xlab("") +
    ylab("")  +
    scale_fill_manual(values = colors) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 15),
          axis.text = element_text(size = 15),
          strip.background =element_rect(fill="white"),
          legend.position = "none")
  return(pp)
}



##ENP
ENPcolor <- c("#1B9E77","#5ab4ac","#c7eae5","#d8b365","#f6e8c3")

##All epochs dadi
ENPSFSdata <- load_dadi_sfs("ENP-44.sfs","ENP"," ENP Data")
ENPSFSmodel <- load_dadi_sfs("ENP.best.sfs","ENP","3 Epoch")
ENPSFSmodel1E <- load_dadi_sfs("ENP.1Epoch.sfs","ENP","1 Epoch")
ENPSFSmodel2E <- load_dadi_sfs("ENP.2Epoch.sfs","ENP","2 Epoch")
ENPSFSmodel4E <- load_dadi_sfs("ENP.4Epoch.sfs","ENP","4 Epoch")

ENP <- rbind(ENPSFSdata,ENPSFSmodel,ENPSFSmodel1E,ENPSFSmodel2E,ENPSFSmodel4E)

##All epochs fsc
ENPSFSdataf <- load_fsc_sfs("ENP_MAFpop0.obs","ENP"," ENP Data",1)
ENPSFSmodelf <- load_fsc_sfs("1D.3Epoch.ENP_MAFpop0.txt","ENP","3 Epoch",2)
ENPSFSmodel1Ef <- load_fsc_sfs("1D.1Epoch_MAFpop0.txt","ENP","1 Epoch",2)
ENPSFSmodel2Ef <- load_fsc_sfs("1D.2Epoch.ENP_MAFpop0.txt","ENP","2 Epoch",2)
ENPSFSmodel4Ef <- load_fsc_sfs("1D.4Epoch.ENP_MAFpop0.txt","ENP","4 Epoch",2)

ENP <- rbind(ENP,ENPSFSdataf,ENPSFSmodelf,ENPSFSmodel1Ef,ENPSFSmodel2Ef,ENPSFSmodel4Ef)
plot_fit_sfs_color_epoch(ENP,"", ENPcolor)



##GOC
GOCcolor <- c("#D95F02","#5ab4ac","#c7eae5","#d8b365","#f6e8c3")

##All epochs dadi
GOCSFSdata <- load_dadi_sfs("GOC-30.sfs","GOC"," GOC Data")
GOCSFSmodel <- load_dadi_sfs("GOC.best.sfs","GOC","2 Epoch")
GOCSFSmodel1E <- load_dadi_sfs("GOC.1Epoch.sfs","GOC","1 Epoch")
GOCSFSmodel2E <- load_dadi_sfs("GOC.3Epoch.sfs","GOC","3 Epoch")
GOCSFSmodel4E <- load_dadi_sfs("GOC.4Epoch.sfs","GOC","4 Epoch")

GOC <- rbind(GOCSFSdata,GOCSFSmodel,GOCSFSmodel1E,GOCSFSmodel2E,GOCSFSmodel4E)



##All epochs fsc
GOCSFSdataf <- load_fsc_sfs("GOC_MAFpop0.obs","GOC"," GOC Data",1)
GOCSFSmodelf <- load_fsc_sfs("1D.3Epoch.GOC_MAFpop0.txt","GOC","3 Epoch",2)
GOCSFSmodel1Ef <- load_fsc_sfs("1D.1EpochGOC_MAFpop0.txt","GOC","1 Epoch",2)
GOCSFSmodel2Ef <- load_fsc_sfs("1D.2Epoch.GOC_MAFpop0.txt","GOC","2 Epoch",2)
GOCSFSmodel4Ef <- load_fsc_sfs("1D.4Epoch.GOC_MAFpop0.txt","GOC","4 Epoch",2)

GOC <- rbind(GOC,GOCSFSdataf,GOCSFSmodelf,GOCSFSmodel1Ef,GOCSFSmodel2Ef,GOCSFSmodel4Ef)
plot_fit_sfs_color_epoch(GOC,"", GOCcolor)
##3*7
pp = ggplot(data = ENP, aes(x = site_freq, y = prop, fill=sfsid)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw()+ facet_grid(vars(sfspop), vars(soft)) +
  xlab("") +
  ylab("")  +
  scale_fill_manual(values = ENPcolor) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.background =element_rect(fill="white"))+
  labs(fill="")
as_ggplot(get_legend(pp))
