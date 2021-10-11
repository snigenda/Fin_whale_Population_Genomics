#! /bin/bash
#$ -wd /u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH
#$ -l h_rt=23:00:00,h_data=30G,highp
#$ -N concat_vcf_oxford
#$ -o /u/project/rwayne/pnunez/reports/oxford.err.txt
#$ -e /u/project/rwayne/pnunez/reports/oxford.out.txt
#$ -m abe
#$ -M pnunez


# This script concat vcf and convert to oxford files
# Author: Meixi Lin, modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub   

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail


#Defining directories ---------------------

workdir=/u/project/rwayne/pnunez/FinWhale/ROHs/RZOOROH

#Main ---------------------------------------

for pop in GOC ENP
do
	Â#New directory
	outdir=$workdir/$pop
  	cd ${outdir}
  	
	#Get all the vcf files made in the previous step
	ls ls *.vcf.gz > vcf_files.txt
	
	#Concat vcf files
 	bcftools concat -f vcf_files.txt -O z -o JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz
	tabix -p vcf JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz
	
	#Convert concat vcf file into Oxford file
	bcftools convert \
	-g JointCalls_08_B_VariantPASS_${pop}_allcontig \
	--chrom --tag GT \
	--threads 8 \
	JointCalls_08_B_VariantPASS_${pop}_allcontig.vcf.gz
  	
done

conda deactivate
