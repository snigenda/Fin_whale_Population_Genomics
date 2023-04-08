#! /bin/bash
#$ -cwd
#$ -l h_rt=00:20:00,h_data=15G
#$ -N easySFSProjection_join
#$ -o /u/project/rwayne/snigenda/finwhale/reports/SFS_projection_names.out.txt
#$ -e /u/project/rwayne/snigenda/finwhale/reports/SFS_projection_names.err.txt
#$ -m abe


# This script concentrate files names of SFS projection per chromosomes
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub SFS_projection_join.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

#Define directories ---------------------------------
workdir=/u/project/rwayne/snigenda/finwhale/SFS/Neutral/SFS_projection_fsc
SFSdir=/u/project/rwayne/snigenda/finwhale/SFS/Neutral/SFS_projection

mkdir -p ${workdir}
mkdir -p ${SFSdir}

#Main ---------------

for pop in GOC ENP
do 
	cd ${workdir}

  	x=$(ls $SFSdir) ; printf "%s\n" "$x" > directories.txt
	while read line; do
		echo -e $SFSdir"/"$line"/fastsimcoal2/"$pop"_MAFpop0.obs" >> $pop"_SFS_projection_files.txt"
	done < directories.txt
	
done 

conda deactivate
