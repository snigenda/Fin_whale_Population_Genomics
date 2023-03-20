#!/bin/bash
#$ -l h_data=20G,h_rt=23:00:00
#$ -wd <homedir>
#$ -o <homedir>/reports/PopStructure/
#$ -e <homedir>/reports/PopStructure/
#$ -m abe

# @version      v2
# @usage        qsub step4_admixture_all50_20210316.sh
# @description  Performs admixture analyses on the LDpruned sites (MAF cutoff = 0.10) on only the ENP individuals
# Author: Paulina Nunez Valencia (pnunez@lcg.unam.mx); Meixi Lin
# Date: Thu Mar 18 15:22:26 2021
# Output:
# 1) 10 runs of admixture from K=1 to K=6
# 2) files with stats
# Notes:
# inputFile may be:
#      - a PLINK .bed file
#      - a PLINK "12" coded .ped file
# 2022-02-08 13:22:00
# Only include ENP individuals

############################################################
## import packages

conda activate gentools

set -eo pipefail

############################################################
## def functions

############################################################
## def variables
DATASET='all50'
REF='Minke'
MAFCUT='10'
TODAY=$(date "+%Y%m%d")

HOMEDIR=<homedir>
WORKDIR=${HOMEDIR}/PopStructure/${DATASET}/${REF}
OUTDIR=${WORKDIR}/revisions_Admixture_ENP_maf10
mkdir -p ${OUTDIR}
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# admixture software source: http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
ADMIXTURE=<software>/admixture/dist/admixture_linux-1.3.0/admixture
# Minke whale reference genome
REFERENCE=<homedir2>/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta

# sample to exclude
EXCLUDE_SAMPLE=("GOC002" "GOC006" "GOC010" "GOC025" "GOC038" "GOC050" "GOC053" "GOC063" "GOC068" "GOC071" "GOC077" "GOC080" "GOC082" "GOC086" "GOC091" "GOC100" "GOC111" "GOC112" "GOC116" "GOC125")
# samples to include (not used)

ENP_SAMPLES="ENPAK19,ENPAK20,ENPAK21,ENPAK22,ENPAK23,ENPAK24,ENPAK25,ENPAK26,ENPAK27,ENPAK28,ENPAK29,ENPAK30,ENPBC16,ENPBC17,ENPBC18,ENPCA01,ENPCA02,ENPCA03,ENPCA04,ENPCA05,ENPCA06,ENPCA07,ENPCA08,ENPCA09,ENPOR10,ENPOR11,ENPOR12,ENPOR13,ENPWA14,ENPWA15"
# Original VCF file
INVCF="${WORKDIR}/JointCalls_all50_filterpass_bialleic_all_LDPruned_maf10_SA_mrF.vcf"
# Out prefix
OUTPREFIX="ENP_LDPruned_maf10_SNPs"

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; GIT commit id ${COMMITID}; Maf cutoff = ${MAFCUT}; Input VCF = ${INVCF}"

############################################################
# start with step3 prune vcf to only include ENP samples
cd $OUTDIR

gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude_sample_name ${ss} "; done) \
--excludeFiltered \
--removeUnusedAlternates \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-V ${INVCF} \
-o ${OUTPREFIX}.vcf

bcftools +counts ${OUTPREFIX}.vcf

############################################################
# convert to plink files
plink --vcf ${OUTPREFIX}.vcf --allow-extra-chr --recode 12 --out ${OUTPREFIX}

############################################################
# start to run admixture
# make a directory to store the P.Q files
mkdir -p Admixture_PQ

echo -e "K,iter,CVERROR,LL" > Admixture_CV_LLsummary_maf${MAFCUT}_${TODAY}.csv
for K in {1..6};do
    for i in {1..10};do
        # -s time the random seed to be generated from the current time
        # --cv In this default setting, the cross-validation procedure will perform 5-fold CV
        ${ADMIXTURE} --cv -s time -j8 ${OUTPREFIX}.ped ${K} | tee log_K${K}.iter${i}.out
        # need to move the files and keep all the runs
        mv ${OUTPREFIX}.${K}.Q Admixture_PQ/${OUTPREFIX}.K${K}.iter${i}.Q
        mv ${OUTPREFIX}.${K}.P Admixture_PQ/${OUTPREFIX}.K${K}.iter${i}.P
        # get the CV error and loglikelihood during each run
        CVERROR=$(awk '/^CV/ {print $4}' log_K${K}.iter${i}.out)
        LL=$(awk '/^Loglikelihood/ {print $2}' log_K${K}.iter${i}.out)
        echo -e "${K},${i},${CVERROR},${LL}" >> Admixture_CV_LLsummary_maf${MAFCUT}_${TODAY}.csv
    done
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
