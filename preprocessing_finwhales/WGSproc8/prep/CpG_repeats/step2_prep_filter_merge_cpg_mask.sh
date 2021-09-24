#!/bin/bash

# Usage: generate cpg island and merged mask
# Wed Mar 18 11:13:03 2020

# This script was executed interactively

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools 

HOMEDIR=/u/project/rwayne/snigenda/finwhale

REFERENCE=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
MASK_CPG=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/ucsc/CpG_ExtUnmasked.bed
MASK_REP=${REFERENCE/.fasta/_repeats.bed}

# use bedtools to merge the two masks 
cd ${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/
cat ${MASK_CPG} > CpG_repeats.bed
cat ${MASK_REP} >> CpG_repeats.bed
mv CpG_repeats.bed CpG_repeats_unmerged.bed
bedtools sort -i CpG_repeats_unmerged.bed > CpG_repeats_sorted.bed
bedtools merge -i CpG_repeats_sorted.bed > CpG_repeats.bed

# Use geneious to check the merge result on one scaffold 
grep "NW_006724458.1" CpG_repeats.bed > ${SCRATCH}/playground/test.bed
cd /u/scratch/m/meixilin/playground
grep "NW_006724458.1" CpG_ExtUnmasked.bed > CpG.bed
grep "NW_006724458.1" GCF_000493695.1_BalAcu1.0_genomic_repeats.bed > repeat.bed

# Fri May 15 10:41:15 2020
# check the rm.out.bed and the CpG_repeats.bed 
RMOUT=GCF_000493695.1_BalAcu1.0_rm.out.bed
MASKREP=GCF_000493695.1_BalAcu1.0_genomic_repeats.bed

# RM out is the repeat masker output file and the MASKREP is the WindowMasker file 
# Nothing showed up because they all failed in the variant filtration 
cd ${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/
head $RMOUT
head $MASKREP
head $REFERENCE

# get the file total length 
awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${RMOUT}
# 969,461,823
awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${MASKREP}
# 741,334,011
awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${MASK_REP}


# get the merge length 
cut -f 1-3 ${RMOUT} > ${SCRATCH}/TwoMasker.out.bed

cat ${MASKREP} >> ${SCRATCH}/TwoMasker.out.bed
bedtools sort -i ${SCRATCH}/TwoMasker.out.bed > ${SCRATCH}/TwoMasker_sorted.out.bed
bedtools merge -i ${SCRATCH}/TwoMasker_sorted.out.bed > ${SCRATCH}/TwoMasker_merged.out.bed
awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${SCRATCH}/TwoMasker_merged.out.bed
# 1,221,985,945

# intersect length 
cut -f 1-3 ${RMOUT} > ${SCRATCH}/RepeatMasker.out.bed
bedtools intersect -a ${SCRATCH}/RepeatMasker.out.bed -b ${MASKREP} > ${SCRATCH}/TwoMasker_intersect.bed
awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${SCRATCH}/TwoMasker_intersect.bed
# 487,884,195

# Final length 
MASK="/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/CpG_repeats.bed"

awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${MASK}
# 768,518,931

MASK2="/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/CpG_repeats_all.bed"

awk -F'\t' '{SUM += $3-$2} END {print SUM}' ${MASK2}
# 1,247,900,490


