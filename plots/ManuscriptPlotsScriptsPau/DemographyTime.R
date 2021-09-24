#Epoch changes
GOC <- data.frame("Time"=c(1000000,176197), "Size"= c(14974,14974))
GOC$pop <- "GOC"

#Epoch changes
ENP <- data.frame("Time"=c(1000000,113985,51), "Size"= c(16481,16481,23921))
ENP$pop <- "ENP"

#bind pops
Pop <- rbind(GOC,ENP)

#Plot log scale timing
ggplot(Pop, aes(x=Size,y=Time, color= pop)) +
  scale_y_log10() +
  annotation_logticks(size = 0.6) + theme_bw() + geom_step(size=2) +
  scale_color_manual(values = c("white","white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line()) +
  ylab("  ") +
  xlab("")

##6*9

ggplot(Pop, aes(x=Size,y=Time, color= pop)) +
  scale_y_log10() +
  annotation_logticks() + theme_bw() + geom_step(size=2) +
  scale_color_manual(values = c("white","white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text=element_text(size=15),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line()) +
  ylab("  ") +
  xlab("")
##5*5