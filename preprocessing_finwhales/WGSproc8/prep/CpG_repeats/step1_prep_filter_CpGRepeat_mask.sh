#! /bin/bash

# Usage: generating bedfiles for variant masking 
# This file was executed interactively 

# 1. The soft masks in NCBI 
cd /u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome
DIR=/u/project/rwayne/jarobins/utils/programs/repeatmask
FASTA=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
flex ${DIR}/code.l && gcc ${DIR}/lex.yy.c && cat ${FASTA} | ./a.out > ${FASTA%.fa*}_repeats.bed
# two files were generated additionally: 
# a.out lex.yy.c, moved them to ./GCF_000493695.1_BalAcu1.0/archive 

# now change the file header names and get rid of the long text annotations: 
# NW_006725630.1 Balaenoptera acutorostrata scammoni unplaced genomic scaffold, BalAcu1.0 scaffold1348, whole genome shotgun sequence --> NW_006725630.1
# here `.*` matches wild card

cd /u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/

sed -i 's/ Balaenoptera acutorostrata.*whole genome shotgun sequence//' GCF_000493695.1_BalAcu1.0_genomic_repeats.bed

# 2. The CpG island masks in UCSC genome browser 
# download the CpG island mask 
# use the unmasked version to combine with the NCBI masks
mkdir -p ucsc
cd ucsc
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/balAcu1/database/cpgIslandExtUnmasked.txt.gz
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/balAcu1/database/chromAlias.txt.gz
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/balAcu1/database/ucscToRefSeq.txt.gz

gunzip chromAlias.txt.gz
gunzip cpgIslandExtUnmasked.txt.gz
gunzip ucscToRefSeq.txt.gz

# modify the cpgIslandExtUnmasked to follow ncbi naming scheme
SCRIPTDIR=/u/project/rwayne/snigenda/finwhale/scripts/WGSproc8/prep/CpG_repeats
WORK=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/ucsc/
Rscript ${SCRIPTDIR}/rename_cpg_ucsc2ncbi.R ${WORK} cpgIslandExtUnmasked.txt ucscToRefSeq.txt CpG_ExtUnmasked.bed


