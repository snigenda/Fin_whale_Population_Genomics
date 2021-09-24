# Title: Format final baleen whale downloads
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jan  8 11:01:31 2021

# preparation --------
rm(list = ls())
cat("\014")
setwd("/Users/linmeixi/google_drive/finwhale/analyses/baleen_genomes/")

library(dplyr)
library(stringr)
# def functions --------
get_SampleName <- function(organism) {
    short = stringr::str_split(organism, pattern = " ")[[1]]
    short = stringr::str_to_title(short)
    short = stringr::str_sub(short, start = 1, end = 3)[1:2]
    short = paste(short, collapse = "")
    return(short)
}

number_SampleName <- function(SampleName) {
    counter = data.frame(unique(SampleName), 0)
    out = c()
    for (ii in 1:length(SampleName)) {
        name = SampleName[ii]
        counter[counter[,1] == name,2] = counter[counter[,1] == name,2] + 1
        out[ii] = paste0(name, str_pad(counter[counter[,1] == name,2], width = 2, pad = "0"), collapse = "")
    }
    return(out)
}
# def variables --------


# load data --------
srarun = read.csv(file = "SraRunTable.txt", stringsAsFactors = FALSE)

# main --------
srainfo = srarun %>%
    dplyr::mutate(Run2 = as.integer(str_remove(Run, "SRR"))) %>% 
    dplyr::arrange(Organism, Run2, ReleaseDate) %>%
    dplyr::mutate(SRA_Accession = Run,
                  R1 = paste0("/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/sra_seq/", Run, "_1.fastq.gz"),
                  R2 = paste0("/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/sra_seq/", Run, "_2.out.fastq.gz"),
                  RGLB = Experiment,
                  RGPU = "NA",
                  RGCN = stringr::str_replace_all(Center.Name, " ", "_"),
                  RGPM = stringr::str_replace_all(Instrument, " ", ""),
                  RGPM = stringr::str_remove(RGPM, "Illumina"),
                  NCBISampleName = Sample.Name,
                  ftpR1 = paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", str_sub(Run,start = 1, end = 6), "/XXX/", Run, "/", Run, "_1.fastq.gz"),
                  ftpR2 = paste0("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/", str_sub(Run,start = 1, end = 6), "/XXX/", Run, "/", Run, "_2.fastq.gz")
                  ) %>%
    dplyr::select(SRA_Accession, R1, R2, starts_with("RG"), Organism, NCBISampleName, AvgSpotLen, Bases, ftpR1, ftpR2)

srainfo$SampleName = number_SampleName(unlist(lapply(srainfo$Organism, get_SampleName)))

srainfo = srainfo %>%
    dplyr::mutate(RGID = stringr::str_c(SampleName, "_A")) %>%
    dplyr::select(SRA_Accession, SampleName, R1, R2, RGID, RGLB, RGPU, RGCN, RGPM, Organism, NCBISampleName, AvgSpotLen, Bases, ftpR1, ftpR2)

write.csv(x = srainfo, file = "/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/baleen_fqpath_fqname.csv", quote = F, row.names = F)

# cleanup --------
