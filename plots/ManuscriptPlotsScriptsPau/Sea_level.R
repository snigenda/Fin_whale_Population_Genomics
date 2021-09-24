#Sea level plot

#Read file (obtain from vaquita paper source)
sea_level <- read.csv("sea_level.csv",header = T)
qq = rep(-100, length(sea_level$RSL_Pmax_m))

  ggplot(data = sea_level, aes(x = time_kyr, y = RSL_Pmax_m)) +
  labs(x = "Years before present", y = "Relative sea level (m)")+
  theme_classic() +
    geom_area(data=data.frame(qq), aes(y=qq), fill='red', alpha=0.5) +
    geom_area(aes(y=RSL_Pmax_m), col='black') 


sea_level$time_kyr <- sea_level$time_kyr*1000  
  
ggplot(sea_level) +
  geom_ribbon(aes(x = time_kyr, ymin = min(RSL_Pmax_m), ymax = RSL_Pmax_m), fill="darkgrey")+
  theme_classic() + scale_x_log10() +
  annotation_logticks()  


ggplot(sea_level) + 
  geom_line(aes(x=time_kyr, y=RSL_Pmax_m), color="black", alpha = 0.5, size = 0.5) +
   scale_x_log10() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())
  

##

