# plot histogram given data
plot_histdt <- function(plotdt, hetpkbmin, hetpkbmax) {
    pph <- ggplot(plotdt, aes(x = hetpkb, color = subpop, fill = subpop)) +
        geom_histogram(breaks = seq(hetpkbmin, hetpkbmax, by = 0.1), closed = 'left') +
        facet_wrap(. ~ sample) 
}


# pull out histgram data from geom_histogram facet variable has to be 'sample'
get_histdt <- function(pph) {
    builddt = ggplot_build(pph)
    # check that the count is the same as ymax
    if (!all(builddt$data[[1]]$count == builddt$data[[1]]$ymax)) {
        stop('Wrong ggplot build')
    }
    # get histogram
    histdt = builddt$data[[1]] %>%
        dplyr::select(PANEL, fill, x, xmin, xmax, count) 
    # get the label of panels
    paneldt = builddt$layout$layout %>%
        dplyr::select(PANEL, sample)
    # merge two dt
    outdt = dplyr::left_join(paneldt, histdt, by = 'PANEL')
    return(outdt)
}

# check this method with previous
otherdt = read.csv(file = '/Users/linmeixi/google_drive/finwhale/analyses/window_het/all50/derive_data/winHet_1Mbwin_1Mbstep_20Per_allsamples_histdt_0.1hetpkbbin_clean_20210622.csv', stringsAsFactors = FALSE)

otherdt2 = reshape2::dcast(data = otherdt, xmax ~ sample, value.var = 'count')
otherdt2 = otherdt2[,-1]
testdt = all50hist01[1:82,3:52]

all(otherdt2 == testdt)
