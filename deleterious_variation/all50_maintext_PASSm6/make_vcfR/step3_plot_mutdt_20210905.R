# Title: FINAL VERSION Plot the master dataset (syn/nonsynTOL/nonsynDEL/LOF) with the invariant sites and bad individuals filtered out
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Feb 18 21:41:45 2021
# Modification: Change layout to adapt for final figure setup 
# Date: Fri Mar 19 18:55:40 2021
# Modification: Change font to ArialMT
# Date: Fri Sep  3 15:35:14 2021



# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

source("~/Lab/fin_whale/scripts_analyses/config/plotting_config.R")

# def functions --------
# https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
scale_override <- function(which, scale) {
    if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
        stop("which must be an integer of length 1")
    }
    
    if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
        stop("scale must be an x or y position scale")
    }
    
    structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
    "CustomFacetWrap", FacetWrap,
    init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
        # make the initial x, y scales list
        scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
        
        if(is.null(params$scale_overrides)) return(scales)
        
        max_scale_x <- length(scales$x)
        max_scale_y <- length(scales$y)
        
        # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
        for(scale_override in params$scale_overrides) {
            which <- scale_override$which
            scale <- scale_override$scale
            
            if("x" %in% scale$aesthetics) {
                if(!is.null(scales$x)) {
                    if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
                    scales$x[[which]] <- scale$clone()
                }
            } else if("y" %in% scale$aesthetics) {
                if(!is.null(scales$y)) {
                    if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
                    scales$y[[which]] <- scale$clone()
                }
            } else {
                stop("Invalid scale")
            }
        }
        
        # return scales
        scales
    }
)

CustomFacetGrid <- ggproto(
    "CustomFacetGrid", FacetGrid,
    init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
        # make the initial x, y scales list
        scales <- ggproto_parent(FacetGrid, self)$init_scales(layout, x_scale, y_scale, params)
        
        if(is.null(params$scale_overrides)) return(scales)
        
        max_scale_x <- length(scales$x)
        max_scale_y <- length(scales$y)
        
        # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
        for(scale_override in params$scale_overrides) {
            which <- scale_override$which
            scale <- scale_override$scale
            
            if("x" %in% scale$aesthetics) {
                if(!is.null(scales$x)) {
                    if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
                    scales$x[[which]] <- scale$clone()
                }
            } else if("y" %in% scale$aesthetics) {
                if(!is.null(scales$y)) {
                    if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
                    scales$y[[which]] <- scale$clone()
                }
            } else {
                stop("Invalid scale")
            }
        }
        
        # return scales
        scales
    }
)


facet_wrap_custom <- function(..., scale_overrides = NULL) {
    # take advantage of the sanitizing that happens in facet_wrap
    facet_super <- facet_wrap(...)
    
    # sanitize scale overrides
    if(inherits(scale_overrides, "scale_override")) {
        scale_overrides <- list(scale_overrides)
    } else if(!is.list(scale_overrides) || 
              !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
        stop("scale_overrides must be a scale_override object or a list of scale_override objects")
    }
    
    facet_super$params$scale_overrides <- scale_overrides
    
    ggproto(NULL, CustomFacetWrap,
            shrink = facet_super$shrink,
            params = facet_super$params
    )
}

facet_grid_custom <- function(..., scale_overrides = NULL) {
    # take advantage of the sanitizing that happens in facet_wrap
    facet_super <- facet_grid(...)
    
    # sanitize scale overrides
    if(inherits(scale_overrides, "scale_override")) {
        scale_overrides <- list(scale_overrides)
    } else if(!is.list(scale_overrides) || 
              !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
        stop("scale_overrides must be a scale_override object or a list of scale_override objects")
    }
    
    facet_super$params$scale_overrides <- scale_overrides
    
    ggproto(NULL, CustomFacetGrid,
            shrink = facet_super$shrink,
            params = facet_super$params
    )
}


# def variables --------
today = format(Sys.Date(), "%Y%m%d")
dataset = 'all50'
ref = 'Minke'
cdstype = 'ALLregions'
prefixlist = c("syn", "nonsynTOL", "nonsynDEL", "LOF")

workdir = paste('/Users/linmeixi/google_drive/finwhale/analyses/DelVar_vcfR/', dataset, ref, sep = '/')
setwd(workdir)

plotdir = './pub_plots/'
dir.create(plotdir)

sink(file = paste0('logs/step3.2_plot_mutdt_', paste(dataset, ref, cdstype, today, sep = '_'), ".log"))
date()
sessionInfo()
# load data --------
# from step4_generate_summarystats_20210218.R output
mutdt = readRDS(file = './derive_data/sum_table/all50_Minke_ALLregions_SUMtable_Seg_PASSm6_Norm_20210219.rds') %>% 
    dplyr::filter(MutType %in% prefixlist) 
mutdt$MutType = factor(mutdt$MutType, levels = prefixlist,
                       labels = c("SYN", "TOL", "DEL", "LOF"))

# main --------
# the allele proportion plot ========
alleledt = reshape2::melt(mutdt, id.vars = c("SampleId", "PopId", "SubPopId", "MutType")) %>%
    dplyr::filter(variable %in% c("normAltAllele"))

pp1 <- ggplot(data = alleledt, aes(x = PopId, y = value, fill = PopId)) +
    geom_boxplot() + 
    scale_fill_manual(values = mycolors) +
    # facet_wrap_custom(. ~ MutType, scales = "free_y",scale_overrides = list(
    #     scale_override(1, scale_y_continuous(breaks = c(25400, 25600, 25800))),
    #     scale_override(2, scale_y_continuous(breaks = c(13700, 13900))),
    #     scale_override(3, scale_y_continuous(breaks = c(5800, 6000))),
    #     scale_override(4, scale_y_continuous(breaks = c(360, 380)))
    # )) +
    facet_wrap(. ~ MutType, scales = "free_y") +
    labs(fill = "Population", x = "Population", y = "Derived Allele Count") +
    theme_pubr() +
    theme(legend.position = "none")

pp1.2 <- pp1 +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox.test",
                       label.x.npc = 'center',
                       vjust = 1,
                       size = 4) +
    theme(text = element_text(family = 'ArialMT'))

ggsave(filename = paste0(plotdir, "MutType_allele_number_box_", today, ".pdf"), plot =  pp1.2, height = 6, width = 6)

# the genotype proportion plot ========
genodt = reshape2::melt(mutdt, id.vars = c("SampleId", "PopId", "SubPopId", "MutType")) %>% 
    dplyr::filter(variable %in% c("normHet", "normHomAlt"))
genodt$variable = as.character(genodt$variable)
genodt$variable = factor(genodt$variable, levels = c("normHomAlt", "normHet"), labels = c("Homozygous Derived", "Heterozygous"))

pp2 <- ggplot(data = genodt,
              aes(x = PopId, y = value, fill = PopId)) +
    geom_boxplot() + 
    scale_fill_manual(values = mycolors) +
    # facet_grid(MutType ~ variable, scale = 'free_y') +
    facet_grid_custom(MutType ~ variable, scale = 'free_y',scale_overrides = list(
        scale_override(1, scale_y_continuous(limits = c(5500, 11500))),
        scale_override(2, scale_y_continuous(limits = c(3000, 6500))),
        scale_override(3, scale_y_continuous(limits = c(1200, 3200)))
        # scale_override(4, scale_y_continuous(breaks = c(100, 150)))
    )) +
    labs(fill = "Population", x = "Population", y = "Genotype Count") +
    theme_pubr() +
    theme(legend.position = "none")

pp2.2 <- pp2 +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox.test",
                       label.x.npc = 'center',
                       vjust = 1,
                       size = 4) +
    theme(text = element_text(family = 'ArialMT'))


ggsave(filename = paste0(plotdir, "MutType_geno_number_box_", today, ".pdf"), plot =  pp2.2, height = 8, width = 5)

# cleanup --------
sink()
closeAllConnections()
