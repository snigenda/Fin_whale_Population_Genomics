#! /bin/bash
#$ -wd /u/home/c/ckyriazi/project-klohmuel/finwhale_sims/output
#$ -l highp,h_rt=250:00:00,h_data=12G
#$ -o /u/home/c/ckyriazi/project-klohmuel/finwhale_sims/output
#$ -e /u/home/c/ckyriazi/project-klohmuel/finwhale_sims/output
#$ -N finwhale_2pop_whalingBott_noMig_071421
#$ -m ae
#$ -t 1:25

# MAKE SURE TO CHANGE NAME OF JOB

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

SLIMDIR=/u/home/c/ckyriazi/project-klohmuel/software/slim_build


${SLIMDIR}/slim /u/home/c/ckyriazi/project-klohmuel/finwhale_sims/scripts/finwhale_2pop_whalingBott_noMig_071421.slim
