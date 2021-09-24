library(ggOceanMaps)


samples <- data.frame("Longitude"=c(-130,-110,-75),"Latitude"=c(46,27,-23),"Pop"=c("ENP","GOC","ESP"))

basemap(limits = c(-140, -67, 58, -58), bathymetry = TRUE,glaciers = T) +
  geom_point(data = samples, aes(x = Longitude, y = Latitude), color="black", size = 6.8) +
  geom_point(data = samples, aes(x = Longitude, y = Latitude, color=Pop), size = 6) +
  theme(legend.position = "none") +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true")+
  scale_color_manual(values=c("#1B9E77","#7570B3","#D95F02")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


             