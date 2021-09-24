#!/bin/bash
#$ -l highp,h_data=8G,h_vmem=INFINITY,h_rt=48:00:00
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/BUSCO_run_20210116.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/BUSCO_run_20210116.err.txt
#$ -m abe

# @version 		v2
# @usage		run the busco analyses
# @description	baleen_genomes/BUSCO_run_20210116
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jan 15 15:09:09 2021
# @modification: Sun Jan 17 11:01:59 2021
# @modification: Run the protein mode first and change the Bryde's rna structure

###########################################################
## import packages
sleep $((RANDOM % 120))

eval "$(/u/project/rwayne/meixilin/miniconda3/bin/conda shell.bash hook)"
conda activate busco

set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
REF=${1} # input the short hand for the genome to compare
BUSCODB=${2} # what busco database to compare it to
PROTEINMODE="TRUE" # default run protein too

TODAY=$(date "+%Y%m%d")
REFDIR=/u/project/rwayne/snigenda/finwhale
HOMEDIR=/u/project/rwayne/meixilin/fin_whale/analyses
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/busco
# WORKDIR=${HOMEDIR}/baleen_genomes/busco
# mkdir -p ${WORKDIR}
mkdir -p ${SCRATCHDIR}

COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)
LOG=${HOMEDIR}/reports/baleen_genomes/BUSCO_${REF}_${BUSCODB:0:5}_${TODAY}.log

###########################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}; The qsub input ${REF} ${BUSCODB}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; git commit id: ${COMMITID}; The qsub input ${REF} ${BUSCODB}" > ${LOG}

cd ${SCRATCHDIR}

# get the transcript and protein files
if [ $REF == 'Minke' ]; then
    TRANSCRIPT=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_rna.fna
    PROTEIN=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_protein.faa
fi
if [ $REF == 'Bryde' ]; then
    TRANSCRIPT=${REFDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni/genome/maker/Balaenoptera_edeni_Balaenoptera_edeni_HiC.fasta_v2.functional.transcripts2N.fasta
    PROTEIN=${REFDIR}/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni/genome/maker/Balaenoptera_edeni_Balaenoptera_edeni_HiC.fasta_v2.functional.proteins.fasta
fi
if [ $REF == 'Fin' ]; then
    TRANSCRIPT=${REFDIR}/cetacean_genomes/fin_whale_genome/GCA_008795845.1_Baphy/GCA_008795845.1_Baphy_rna_from_genomic.fna
    PROTEIN=${REFDIR}/cetacean_genomes/fin_whale_genome/GCA_008795845.1_Baphy/GCA_008795845.1_Baphy_protein.faa
fi
if [ $REF == 'Blue' ]; then
    TRANSCRIPT=${REFDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_rna.fna
    PROTEIN=${REFDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_protein.faa
fi
if [ $REF == 'Humpback' ]; then
    TRANSCRIPT=${REFDIR}/cetacean_genomes/humpback_whale_genome/GCA_004329385.1_megNov1/megNov1_transcripts.final.fasta
    PROTEIN=${REFDIR}/cetacean_genomes/humpback_whale_genome/GCA_004329385.1_megNov1/megNov1_proteins.final.fasta
fi

# the two other rna files
if [ $REF == 'Minke2' ]; then
    PROTEINMODE="FALSE"
    TRANSCRIPT=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_rna_from_genomic.fna
fi
if [ $REF == 'Blue2' ]; then
    PROTEINMODE="FALSE"
    TRANSCRIPT=${REFDIR}/cetacean_genomes/blue_whale_genome/GCF_009873245.2_mBalMus1.pri.v3/GCF_009873245.2_mBalMus1.pri.v3_rna_from_genomic.fna
fi

# should you run the protein mode?
if [ ${PROTEINMODE} == "TRUE" ]
then
    # run the protein modes (note change from default setting, e-value cutoff)
    if [ -d busco_downloads/lineages/${BUSCODB} ]
    then
        busco -i ${PROTEIN} \
        -l ${BUSCODB} \
        --offline \
        -c 8 \
        -e 1e-6 \
        -m prot \
        -o ${REF}_prot_${BUSCODB:0:5} &>> ${LOG}
    else
        busco -i ${PROTEIN} \
        -l ${BUSCODB} \
        -c 8 \
        -e 1e-6 \
        -m prot \
        -o ${REF}_prot_${BUSCODB:0:5} &>> ${LOG}
    fi
fi

# run the transcript modes (note change from default setting, e-value cutoff)
if [ -d busco_downloads/lineages/${BUSCODB} ]
then
    busco -i ${TRANSCRIPT} \
    -l ${BUSCODB} \
    --offline \
    -c 8 \
    -e 1e-6 \
    -m tran \
    -o ${REF}_tran_${BUSCODB:0:5} &>> ${LOG}
else
    busco -i ${TRANSCRIPT} \
    -l ${BUSCODB} \
    -c 8 \
    -e 1e-6 \
    -m tran \
    -o ${REF}_tran_${BUSCODB:0:5} &>> ${LOG}
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}


###########################################################
# NOTES and tests

# For memory issues, see this:
# https://gitlab.com/ezlab/busco/-/issues/356
# Genome ftp: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/

# BUSCODB="cetartiodactyla_odb10"
# REF="Minke"
# REFDIR=/u/project/rwayne/snigenda/finwhale
# TRANSCRIPT=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_rna.fasta
# PROTEIN=${REFDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_protein.faa
# LOG="test.log"

# # functions to gunzip and keep the original gz file
# cpandgunzip() {
# local FILE=${1}
# cp ${FILE} ${FILE/.gz/2.gz}
# gunzip ${FILE/.gz/2.gz}
# mv ${FILE/.gz/2} ${FILE/.gz}
# }

