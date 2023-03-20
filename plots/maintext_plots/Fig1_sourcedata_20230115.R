# Title: Generate source data for Fig 1a, the maps
# Author: Meixi Lin
# Date: Sun Jan 15 20:04:52 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("<homedir>/fin_whale/FinWhale_PopGenomics_2021")

library(readxl)
library(dplyr)
# def functions --------

# def variables --------

# load data --------
enpdt = read_excel(path = '~/Google Drive/My Drive/finwhale/preprocessing/sample_info/ENP_whales_samples_INFO_plot.xlsx')[-31,]
gocdt = read_excel(path = '~/Google Drive/My Drive/finwhale/preprocessing/sample_info/Metadata_GOC_samples.xlsx')

# main --------
enpdt = enpdt %>%
    dplyr::rename(sampleid = SEQUENCING_ID, location = COAST, longitude = Longitude, latitude = Latitude) %>%
    dplyr::mutate(pop = 'ENP') %>%
    dplyr::select(sampleid, location, longitude, latitude, pop)

gocdt = gocdt %>%
    dplyr::mutate(pop = 'GOC',
                  sampleid = stringr::str_remove(WGS_ID,'-')) %>%
    dplyr::select(sampleid, location, longitude, latitude, pop) %>%
    dplyr::mutate(longitude = ifelse(longitude > 0, -longitude, longitude))

outdt = rbind(enpdt, gocdt) %>%
    dplyr::arrange(sampleid)

# redact latlong
outdt2 = outdt %>%
    dplyr::mutate(longitude = 'Available upon request',
                  latitude = 'Available upon request')

# output files --------
write.csv(outdt, file = './source_data/Fig1a-wLatLon.csv')
write.csv(outdt2, file = './source_data/Fig1a.csv')

# cleanup --------
date()
closeAllConnections()
