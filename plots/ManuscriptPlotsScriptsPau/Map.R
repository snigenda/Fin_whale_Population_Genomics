library(ggOceanMaps)
library(ggsn)
library(shades)
library(ggpubr)
library(ggforce)
library(ggmap)
library(cowplot)


setwd("C:/Users/gabri/Documents/MorenoLab/Ballenas/Minke_all48/Map")


#Data 
ENP <- read.csv("ENP_whales_samples_INFO.csv")
ENP <- ENP[,c(7,8,13)]


ENP_loc <- basemap(data=ENP,bathymetry = TRUE)  + 
  geom_point(data = ENP, aes(x = Longitude, y = Latitude, color=pop), size=8, color="#148563") +
  theme(legend.position = "none") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", which_north = "true")+
  theme(text = element_text(size=25))

GOC <- read.csv("Metadata_GOC_samples.csv")
GOC <- GOC[,c(3,4)]
GOC <- GOC[-c(16,19),]
colnames(GOC) <- c("Latitude", "Longitude")
GOC$COAST <- "GOC"

GOC_loc <- basemap(data = GOC, bathymetry = TRUE, limits = c(-118, -107, 33, 22))  + 
  geom_point(data = GOC, aes(x = Longitude, y = Latitude, color=pop), size=8, color="#D95F02") +
  geom_point(data = GOC, aes(x = Longitude, y = Latitude), size=8, color="black", pch=21, stroke=2) +
  theme(legend.position = "none") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", which_north = "true")+
  theme(text = element_text(size=25))

samples <- rbind(ENP,GOC)
colnames(samples) <- c("lat","lon","pop")

all_pop <- basemap(data = samples, bathymetry = TRUE, glaciers = TRUE)  + 
  geom_point(data = samples, aes(x = lon, y = lat, color=pop), size= 4) +
  theme(legend.position = "none") +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true")+
  scale_color_manual(values=loccolors)

###########



loccolors = c("#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")


#Data 
ENP <- read.csv("ENP_whales_samples_INFO.csv")
ENP <- ENP[,c(7,8,13,15)]
GOC <- read.csv("Metadata_GOC_samples.csv")
GOC <- GOC[,c(3,4)]
GOC <- GOC[-c(16,19),]
colnames(GOC) <- c("Latitude", "Longitude")
GOC$COAST <- "GOC"

ENP_loc <- basemap(data=ENP,bathymetry = TRUE, grid.col = NA, 
                   limits = c(-174, -103, 0, 65), land.col = "#BBBABA")  + 
  geom_point(data = ENP, aes(x = Longitude, y = Latitude, color=COAST), size=6) +
  theme(legend.position = "none") +
  annotation_scale(location = "br") +
  annotation_north_arrow(location = "tr", which_north = "true")+
  theme(text = element_text(size=15)) +
  scale_colour_manual(values = loccolors) +
  labs(x = NULL, y = NULL)
  

GOC_loc <- basemap(data = GOC, bathymetry = TRUE, grid.col = NA,
                   limits = c(-118, -108, 33, 22), land.col = "#BBBABA")  + 
  geom_point(data = GOC, aes(x = Longitude, y = Latitude, color=COAST), size=3.5, color="#D95F02")  +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x = NULL, y = NULL)


ggdraw(ENP_loc) + draw_plot(GOC_loc, x = 0.16, y = 0.06, width = 0.42, 
                            height = 0.42) 

## 5*5

##5*5