library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)

dir <- "C:/Users/gabri/Documents/MorenoLab/Ballenas/Results/Demography"

setwd(dir)


#Data
sfs <- scan("ENP-GOC.sfs", skip = c(1,3))
sfs <- matrix(sfs, ncol = 45, nrow = 31)
sfs[1,1] <- NA
sfs2 <- melt(sfs)
sfs2$value[sfs2$value==0] <- NA
sfs2$Type <- "Data"


#Model
model <- scan("ENP-GOC.best.sfs", skip = c(1,3))
model <- matrix(model, ncol = 45, nrow = 31)
model <- model * 1000000
model[1,1] <- NA
model2 <- melt(model)
model2$value[model2$value==0] <- NA
model2$Type <- "Model"


#Join

forplot <- rbind(sfs2, model2)

#Baseplot
heatmap(sfs2,Colv = NA, Rowv = NA)


##Gpplot colors

ggplot(forplot, aes(x=Var1,y=Var2,fill=value)) +
  geom_tile() + facet_wrap(~Type)+
  scale_fill_viridis(discrete=FALSE, direction = -1, na.value="white", trans = "log10")  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.position = "none",
        strip.text.x = element_text(size = 15)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_fixed() +
  xlab("GOC") + ylab("ENP")

##4.5*6

##Gpplot gray

ggplot(forplot, aes(x=Var1,y=Var2,fill=value)) +
  geom_tile() + facet_wrap(~Type)+
  scale_fill_viridis(discrete=FALSE, direction = -1, na.value="white", trans = "log10",option = "W")  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.position = "none",
        strip.text.x = element_text(size = 15)) +
  coord_fixed() +
  xlab("GOC") + ylab("ENP")



##Lege3d

pp <- ggplot(forplot, aes(x=Var1,y=Var2,fill=value)) +
  geom_tile() + facet_wrap(~Type)+
  scale_fill_viridis(discrete=FALSE, direction = -1, na.value="white", trans = "log10")  +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 15),
        legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 7)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_fixed() +
  xlab("GOC") + ylab("ENP") +
  labs(fill="")

as_ggplot(get_legend(pp))

