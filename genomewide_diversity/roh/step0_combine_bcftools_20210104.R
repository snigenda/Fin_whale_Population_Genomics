# Title: Combine bcftools output
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Jan  6 16:59:51 2021

# preparation --------
rm(list = ls())
cat("\014")

# def functions --------
read_RG_roh <- function(data){
    output <- read.table(file = data, col.names=c("row_type","sample","chrom","start","end","length","num_markers","qual"), fill=T, stringsAsFactors = FALSE)
    output1 <- base::subset(output, row_type == "RG")
    return(output1)
}

# def variables --------
outdir = "/u/project/rwayne/meixilin/fin_whale/analyses/ROH/rohBcftools/"
today = format(Sys.Date(), "%Y%m%d")

# load data --------
ind50 = read.csv(file = "/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/popmap_all50.csv", stringsAsFactors = F, header = 1)

# main --------
# read files
rglines = data.frame()
for (ii in 1:nrow(ind50)) {
    sample=ind50[ii,'SampleId']
    pop=ind50[ii,'PopId']
    filename = paste0("/u/project/rwayne/snigenda/finwhale/analyses_results/structure_diversity/roh/rohBcftools/",
                      pop, "/",sample, "_", pop, "_concat_fwhale_roh_bcftools_G30_ACANGT", ".out.gz")
    if(!file.exists(filename)) {
        print(paste("WARNING:", sample, "No file named", filename))
    } else {
        temp = read_RG_roh(filename)
        rglines = base::rbind(rglines, temp, stringsAsFactors = FALSE)
    }
}

# output lines
write.csv(rglines, file = paste0(outdir, "all50m3_concat_fwhale_roh_bcftools_G30_ACANGT_", today, ".csv"))

# cleanup --------
closeAllConnections()