#! /bin/bash

# Usage: snpEff database building only for the Minke whale 
# Less error compared with before using gtf files 

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

REF="Minke"
Database=Baac01.10776 # 10776 is the scaffold number
Snpeffdir=/u/project/rwayne/software/finwhale/miniconda2/envs/gentools/share/snpeff-4.3.1t-2
# Snpeffdir=/u/home/m/meixilin/snpEff/snpEff
HOMEDIR=/u/project/rwayne/snigenda/finwhale
MINKEDIR=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0

# chooses the reference genome and gtf file to construct the SnpEff database
REFERENCE=${MINKEDIR}/GCF_000493695.1_BalAcu1.0_genomic.fasta
GTF=${MINKEDIR}/GCF_000493695.1_BalAcu1.0_genomic.gtf

# the config files already modified manually 
# vi snpEff.config 
# append to file 
# % Baac01.10776.genome : Balaenoptera_acutorostrata_scammoni
# % 	# Baac01.10776.reference : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/493/695/GCF_000493695.1_BalAcu1.0/

# goes to the data directory within the SnpEff directory and puts the necessary files to create the SnpEff database in place with the appropriate names
mkdir -p ${Snpeffdir}/data/${Database}
mkdir -p ${Snpeffdir}/data/genomes
cd ${Snpeffdir}/data/${Database}

ln -s ${GTF} genes.gtf

cd ${Snpeffdir}/data/genomes
ln -s ${REFERENCE} ${Database}.fa


# NOW STARTS BUILDING ########
# Use auto detection for file format 
cd ${Snpeffdir}
java -Xmx10g -jar snpEff.jar build -debug ${Database} &> ${Database}.build.log
##############################
