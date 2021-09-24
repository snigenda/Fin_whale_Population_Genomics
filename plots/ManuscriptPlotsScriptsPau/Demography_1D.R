#plot 1D demographic models as lines

GOC <- data.frame("Time"=c(1000000,176197,0), "Size"= c(14974,14974,5304))
GOC$pop <- "GOC"

ENP <- data.frame("Time"=c(1000000,113985,51,0), "Size"= c(16481,16481,23921,306))
ENP$pop <- "ENP"

Pop <- rbind(GOC,ENP)

ggplot(Pop, aes(x=Time,y=Size, color= pop)) +
  scale_x_log10() +
  annotation_logticks() + theme_light() + geom_step(size=2) +
  scale_color_manual(values = c("#1B9E77","#D95F02"))

