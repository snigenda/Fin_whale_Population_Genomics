#! /bin/bash
#$ -wd /u/project/rwayne/snigenda/finwhale/scripts/get_neutral_regions
#$ -l h_rt=04:00:00,h_data=5G,h_vmem=10G
#$ -o /u/project/rwayne/snigenda/finwhale/reports/exon_dist.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/exon_dist.err.txt
#$ -m abe

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

HOMEDIR=/u/project/rwayne/snigenda/finwhale
WORKDIR=${HOMEDIR}/Neutral_regions
REFDIR=${HOMEDIR}/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0
GFF=${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic.gff
CDS_REGIONS=${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic_cds.bed
EXONS=${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic_exons.bed
CDRS=${REFDIR}/GCF_000493695.1_BalAcu1.0_allCodingRegions_merged.bed
hqSites=${WORKDIR}/HiQualCoords/all_HQCoords_sorted_merged.bed

mkdir -p ${WORKDIR}/DistanceFromExons

### Get exonic regions or coding regions (do only once) and make sure is sorted

# awk '$3=="exon"' ${GFF} | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > ${REFDIR}/GCF_000493695.1_BalAcu1.0_genomic_exons.bed
## There are discrepancies between the two methods, Annabel's approach results in 499895 lines, while mine has 499380 lines.
## I ended up using mine because it seems more specific to the column were the features are defined in gff files.

bedtools closest -d -a ${hqSites} -b ${EXONS} > ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt
bedtools closest -d -a ${hqSites} -b ${CDRS} > ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt


# Exploring and choosing different distances

# The last column (8) is the distance; I want it to be at least 10,000, and want to keep track of the distance. Collect all that are >10,000 away.
# Pick the ones with specific distance (awk) from exons, we tried 3 different distances: 10Kb, 20Kb, 50Kb.
awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromExons.0based.bed
awk -F'\t' '{OFS="\t";if($8>20000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromExons.0based.bed
awk -F'\t' '{OFS="\t";if($8>50000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min50kb_DistFromExons.0based.bed
# Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)

awk -F'\t' '{OFS="\t";if($7>10000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromCDR.0based.bed
awk -F'\t' '{OFS="\t";if($7>20000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min20kb_DistFromCDR.0based.bed
awk -F'\t' '{OFS="\t";if($7>50000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min50kb_DistFromCDR.0based.bed

> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
for i in `ls ${WORKDIR}/DistanceFromExons/*DistFromExons.0based.bed`
do
echo $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
done

> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
for i in `ls ${WORKDIR}/DistanceFromExons/*DistFromCDR.0based.bed`
do
echo $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
done

conda deactivate

