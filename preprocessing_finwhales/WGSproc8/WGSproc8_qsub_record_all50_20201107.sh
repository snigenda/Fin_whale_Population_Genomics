###########################################################
# Sat Nov  7 23:10:24 2020
# Get the statistics again
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc8
qsub -t 96 prep/vcf_stats/prep_filter_vcf_getdp_fast.sh meixilin Minke

# Sun Nov  8 13:18:15 2020
# check if the statistics is the same as before. YES!
# check if the reports had no error. YES!
qsub -t 1-95 prep/vcf_stats/prep_filter_vcf_getdp_fast.sh meixilin Minke


###########################################################
# Mon Nov  9 23:15:53 2020
# Perform the first set of filters
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc8
qsub -t 1-96 WGSproc8_a_FilterVCFfile_GATKfilter_20201109.sh meixilin Minke

# check if all jobs completed without error
qacct -j 5217468 | grep 'exit_status' | grep -c 'exit_status  0' # all 96 jobs finished without error

cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke/logs
for ii in {01..96}
do
if ! grep "Done WGSproc8_a_GATKfilter for all50 Minke" WGSproc8_Minke_MarkDuplicates_${ii}_progress_all50.log; then
echo "$ii is not done"
fi
done

###########################################################
# # Tue Nov 10 09:00:11 2020
# Perform the second set of filters
cd /u/project/rwayne/snigenda/finwhale/scripts/WGSproc8
qsub -t 1-96 WGSproc8_b_FilterVCFfile_customfilter_20201109.sh meixilin Minke

# Tue Nov 10 19:26:52 2020
# check if all jobs completed without error
qacct -j 5228838 | grep 'exit_status' | grep -c 'exit_status  0'

cd /u/project/rwayne/snigenda/finwhale/filteredvcf/all50/Minke/logs
for ii in {01..96}
do
if ! grep "Done WGSproc8_b_customfilter" WGSproc8_Minke_MarkDuplicates_${ii}_progress_all50.log; then
echo "$ii is not done"
fi
done

# DONE!