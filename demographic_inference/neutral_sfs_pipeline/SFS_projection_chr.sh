#! /bin/bash
#$ -cwd
#$ -l h_rt=05:00:00,h_data=5G
#$ -N easySFSProjetion_chr
#$ -o /u/project/rwayne/snigenda/finwhale/reports/SFS_projection.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/SFS_projection.err.txt
#$ -t 1-96
#$ -m abe


# This script runs SFS projection  per chromosomes
# Author: Meixi Lin, modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub SFS_projection_chr.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

#Define directories ---------------------------------
homedir=/u/project/rwayne/snigenda/finwhale
workdir=${homedir}/SFS/Neutral
scriptdir=${homedir}/scripts/SFS
easySFS=${scriptdir}/easySFS_a.py 
vcfdir=${workdir}/SNPs
IDX=$(printf %02d ${SGE_TASK_ID})


#Main ---------------

cd ${workdir}
popfile=${homedir}/scripts/config/pop_map.txt
outdir=${workdir}/SFS_projection

mkdir -p $outdir
cd ${outdir}

python $easySFS -i ${vcfdir}/SNPs_neutral_for_SFS_${IDX}.vcf.gz -p ${popfile} --proj 44,30 -a -f -v -o SNPs_easySFS_projection_${IDX} -maxHetFilter 0.75 

conda deactivate

