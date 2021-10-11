#! /bin/bash
#$ -wd /u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH
#$ -l h_rt=23:00:00,h_data=20G,highp
#$ -N filter_vcftools
#$ -o /u/project/rwayne/pnunez/reports/filter_vcf_test.err.txt
#$ -e /u/project/rwayne/pnunez/reports/filter_vcf_test.out.txt
#$ -t 1-96
#$ -M pnunez


# This script filters vcf and subset them per pops
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub   

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail


#Defining directories ---------------------

workdir=/u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH
vcfdir=/u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke
IDX=$(printf %02d ${SGE_TASK_ID})


#Main ---------------------------------------

for pop in GOC ENP
do
  outdir=$workdir/$pop
  mkdir -p $outdir
  cd ${outdir}

  vcftools --gzvcf ${vcfdir}/JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz \
  --remove-filtered-all \
  --keep ${pop}_samples.txt \
  --recode \
  --recode-INFO-all \
  --out ${workdir}/${pop}/VariantFiltration_${IDX}_OnlyPass.${pop}

  gzip ${workdir}/${pop}/VariantFiltration_${IDX}_OnlyPass.${pop}.recode.vcf

done

conda deactivate
