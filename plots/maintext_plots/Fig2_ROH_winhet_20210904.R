# Title: FINAL Fig 2. ROH and window het plots
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Sep  6 22:05:10 2021

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd('/Users/linmeixi/Lab/finwhale_manuscript/')

library(dplyr)
library(ggplot2)
library(ggpubr)
# def functions --------
get_order <- function(genomehet, rohsum, rohlevels) {
    sampleorder = genomehet %>% dplyr::arrange(GenomeHet)
    sampleorder = sampleorder$SampleId
    genomehet$SampleId = factor(genomehet$SampleId, levels = sampleorder)
    rohsum$sample = factor(rohsum$sample, levels = sampleorder)
    rohsum$rohcat = factor(rohsum$rohcat, levels = rohlevels)
    return(list(genomehet, rohsum))
}

plot_indv_winhet <- function(winhetdt, samplename, indivslab) {
    het <- winhetdt[winhetdt$sample == samplename,]
    pphet <- ggplot() +
        geom_col(data = het, mapping = aes(x = merge_start, y = hetpkb, color = subpop, fill = subpop), size = 1e-6) +
        scale_color_manual(values = loccolors) +
        scale_fill_manual(values = loccolors) +
        annotate("text", x = 3.8e+08, y = 7.2 , 
                 label = indivslab[samplename], size = 7.5, family = 'ArialMT') +
        facet_wrap(. ~ sample, nrow = 1, ncol = 1) +
        coord_cartesian(ylim=c(0,7.8), expand = FALSE) + 
        scale_x_continuous(breaks = c(0,0.5,1,1.5,2)*10^9,
                           labels = c('0', '5.0e+8','1.0e+9', '1.5e+9', '2.0e+9')) +
        theme_bw() +
        theme(legend.position = 'none', 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text=element_text(size=20),
              strip.text.x = element_text(size = 25),
              strip.background = element_rect(fill="white")) +
        theme(text = element_text(family = 'ArialMT'))
    
    # used a different than pphet axis lim (7.8 --> 7) (still included all data but gives more space)
    pphist <- ggplot(het, aes(x = hetpkb, color = subpop, fill = subpop)) +
        geom_histogram(binwidth = 0.15,size=0.1) +
        facet_wrap(. ~ sample, nrow = 1, ncol = 1) +
        scale_color_manual(values = loccolors) +
        scale_fill_manual(values = loccolors) +
        scale_y_continuous(expand = c(0,1)) +
        scale_x_continuous(breaks = c(0,2,4,6)) +
        coord_cartesian(xlim = c(0,7), ylim = c(0,450)) +
        theme_bw() +
        theme(legend.position = 'none',
              strip.text = element_blank(),
              plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text=element_text(size=20)) +
        theme(axis.line = element_line(color = 'black')) +
        theme(text = element_text(family = 'ArialMT'))
    
    outpp <- ggarrange(pphet, NULL, pphist, nrow = 1, widths =  c(4,0.12,1))
    return(outpp)
}

get_multileg <- function(pp) {
    # get the color scale
    pp1 <- pp +
        guides(fill = 'none')
    colleg <- as_ggplot(get_legend(pp1))
    # get the fill scale
    pp2 <- pp +
        guides(color = 'none')
    fillleg <- as_ggplot(get_legend(pp2)) 
    return(list(colleg, fillleg))
}

# def variables --------
plotdir = './plots/maintext/'
dir.create(plotdir)
# roh colors and categories ========
nrohcat = 3
rohcatbrewer = c("#B8D8D8","#7A9E9F","#4F6367")
rohlens = c(0.1,1,5)*1e+6
rohcats = c('0.1_1','1_5','5_Inf')
lenslab <- c('[0.1, 1) Mb', '[1, 5) Mb', '[5, Inf) Mb')
names(rohcatbrewer) = rohcats

source('./scripts/config/plotting_config.R')

# load data --------
zooroh_plotdt = readRDS(file = './data/roh/derived_data/rohsummarylong_zoobcf_20210908.rds') %>%
    dplyr::filter(software == 'zoo')
genomehet = read.csv(file = "./data/genome_stats/derived_data/all50_genomewide_heterozygosity_20210824.csv",
                     row.names = 1, stringsAsFactors = F)
othergenomehet = read.csv(file = "./scripts/config/Baleen_Genomewide_Het_20210906.csv", stringsAsFactors = F)

# windowed het ========
winhetdt = readRDS(file = './data/window_het/derived_data/winHet_1Mbwin_1Mbstep_20Per_all50_plotdt_20210824.rds')

indivs = c('ENPAK19','ENPCA09','GOC002','GOC125')
indivslab = c("Mean het. = 1.758","Mean het. = 1.773","Mean het. = 1.035", "Mean het. = 1.195")
names(indivslab) = indivs
winhetdt = winhetdt %>% 
    dplyr::filter(sample %in% indivs)

# # genomehet for plot (how to derive indivslab)
# genomehet_forplot <- genomehet %>%
#     dplyr::filter(SampleId %in% indivs) %>%
#     dplyr::mutate(GenomeHetText = paste('Mean het. =', sprintf("%.3f", round(1000*GenomeHet, digits = 3)))) %>%
#     dplyr::rename(sample = SampleId)

# # get axis limit
# max(winhetdt$hetpkb) # 6.76

# main --------
# preformat the zooroh_plotdt
zooroh_plotdt = zooroh_plotdt %>%
    dplyr::mutate(totalROHMB = sumcat/1e+6)
# sort by genomehet
forplot=get_order(genomehet, zooroh_plotdt, rohlevels = rev(rohcats))

# plotting ROH ========
p1 <-ggplot() +
    scale_y_continuous(limits = c(0, 2.6), name = "", expand= c(0,0),
                       sec.axis = sec_axis(~ . * 5e+2, name="  ")) +
    geom_bar(data = forplot[[2]], 
             aes(x = sample, y = totalROHMB/5e+2, fill = factor(rohcat)),
             stat="identity",position="stack") + 
    geom_point(data = forplot[[1]], aes(x = SampleId, y = GenomeHet*1e+3, color = PopId), size = 2.5)+
    geom_point(data = forplot[[1]], aes(x = SampleId, y = GenomeHet*1e+3),colour = "grey90", size = 1)+
    geom_point(data = othergenomehet,x=-0.7, aes(y = Observed.pi*1e+3), color = "black", shape = 1, size =2.5)+
    geom_text(data = othergenomehet, x = -1.5, aes(y = Observed.pi*1e+3, label = Common.Name), hjust = 'right', color = "black", angle = 45, size= 3)+
    scale_fill_manual(values = rohcatbrewer, labels = rev(lenslab)) +
    coord_flip(xlim = c(0,50), clip = 'off') +
    scale_color_manual(values = mycolors)+
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x.top = element_text(size=10),
          axis.title.x.bottom = element_text(size=10),
          axis.text.y = element_text(size=11),
          axis.text.x.top = element_text(size=12),
          axis.text.x.bottom = element_text(size=12),
          # legend.position = 'none',
          legend.title = element_blank(),
          plot.margin = unit(c(-.1,1,3,1), "lines")) +
    annotate("text", x=-8.5,y=1.1,label="  ",size= 3.5, family = 'ArialMT') +
    theme(text = element_text(family = 'ArialMT'))

# get legend
rohleg <- get_multileg(p1)

# remove legend
p1.1 <- p1 +
    theme(legend.position = 'none')

# plotting window het ========
p2list <- lapply(indivs, plot_indv_winhet, winhetdt = winhetdt, indivslab = indivslab)

# output files --------
ggsave(filename = 'Fig2A_hetlegend.pdf', plot = rohleg[[1]], path = plotdir, height = 2, width = 2)
ggsave(filename = 'Fig2A_rohlegend.pdf', plot = rohleg[[2]], path = plotdir, height = 2, width = 2)
ggsave(filename = 'Fig2A_content.pdf', plot = p1.1, path = plotdir, height = 9, width = 5.2)

for (ii in 1:length(indivs)) {
    ggsave(filename = paste0('Fig2B_content_', indivs[ii],'.pdf'), 
           plot = p2list[[ii]], path = plotdir, height = 4, width = 10.5)
}


# cleanup --------
date()
closeAllConnections()
