#! /bin/bash
#$ -cwd
#$ -l h_rt=23:00:00,h_data=15G
#$ -N easySFSProjetion_monomorphic
#$ -e /u/project/rwayne/snigenda/finwhale/reports/SFS_monomorphic.err.txt
#$ -o /u/project/rwayne/snigenda/finwhale/reports/SFS_monomorphic.out.txt
#$ -t 1-96
#$ -m abe


# This script calculates monomorphic site.
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub SFS_add_monomorphic_sites.sh

source /u/local/Modules/default/init/modules.sh
module load python/2.7

set -o pipefail

#Define directories ---------------------------------

homedir=/u/project/rwayne/snigenda/finwhale
workdir=${homedir}/SFS/Neutral
vcfdir=${homedir}/Neutral_regions/neutralVCFs
outdir=${workdir}/SFS_projection_Monomorphic
monoscript=${homedir}/scripts/SFS/getMonomorphicProjectionCounts.1D.2DSFS.py
IDX=$(printf %02d ${SGE_TASK_ID})


#Main ---------------

popfile=${homedir}/scripts/config/pop_map.txt
projections="44,30" 

mkdir -p ${outdir}
cd ${outdir}

########## get counts of monomorphic sites to add to the SFSes ############

python ${monoscript} --vcf ${vcfdir}/Neutral_sites_SFS_${IDX}.vcf.gz --popMap ${popfile} --proj ${projections} --popIDs ENP,GOC --outdir ${outdir} --outPREFIX ${IDX}

