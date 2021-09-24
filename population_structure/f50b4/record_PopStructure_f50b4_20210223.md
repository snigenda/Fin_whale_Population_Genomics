# step0: generate SNPRelate files

```bash
# Tue Feb 23 22:22:25 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub PopStructure/f50b4/wrapper_step0_vcf2gds_f50b4_20210223.sh

qacct -j 6651232

# Wed Feb 24 22:39:26 2021
# Fix the name later
cd /u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke/

mv -v JointCalls_f50b4_08_B_VariantFiltration_all.gds JointCalls_f50b4_08_B_VariantFiltration_bialleic_all.gds

echo -e "[$(date "+%Y-%m-%d %T")] Rename File JointCalls_f50b4_08_B_VariantFiltration_all.gds to JointCalls_f50b4_08_B_VariantFiltration_bialleic_all.gds" >> /u/project/rwayne/meixilin/fin_whale/analyses/reports/PopStructure/step0_vcf2gds_f50b4_20210223.out.txt
```

# step1: gds LD pruning

```bash
# Thu Feb 25 16:29:06 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

qsub PopStructure/f50b4/wrapper_step1_gdsLDPruning_f50b4_20210223.sh

qacct -j 6667127

cd /u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke
gzip JointCalls_f50b4_filterpass_bialleic_all_20210226.txt
```
# step1.x: transfer the files locally

```bash
# Fri Feb 26 17:28:10 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke/

rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" --exclude 'JointCalls_f50b4_08_B_VariantFiltration_bialleic_all.gds' ${REMOTE} ${LOCAL}

REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke/JointCalls_f50b4_08_B_VariantFiltration_bialleic_all.gds
rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}

```

# step2: gds file to vcf and plink formats

```bash
```

# step 3.1: run APE and bootstrap (not plottings)

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
qsub -N 'ape_mafNA' PopStructure/f50b4/wrapper_step3.1_ApePhylogeny_f50b4_hoff_20210301.sh 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA.gds' 'f50b4_pass_bialleic_all_LDPruned_mafNA'
qsub -N 'ape_maf05' PopStructure/f50b4/wrapper_step3.1_ApePhylogeny_f50b4_hoff_20210301.sh 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf05.gds' 'f50b4_pass_bialleic_all_LDPruned_maf05'
qsub -N 'ape_maf10' PopStructure/f50b4/wrapper_step3.1_ApePhylogeny_f50b4_hoff_20210301.sh 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf10.gds' 'f50b4_pass_bialleic_all_LDPruned_maf10'

# JOB ID: 6718915-6718917
qacct -j 6718917 | grep 'exit_status' | grep -c 'exit_status  0'

cd /u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke

```

# step 3.2: plot the APE trees

```bash
# Tue Mar  2 09:49:16 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/PopStructure/f50b4/Minke/

rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" --delete --update ${REMOTE} ${LOCAL}

rsync -ahv -e "ssh -i ~/.ssh/id_hoff_rsa" --update ${REMOTE} ${LOCAL}


# run the plotting
# Wed Mar 10 13:33:28 2021
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/f50b4/step3.2_ApePhylogeny_f50b4_plot_20210301.R '/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke' 'f50b4_pass_bialleic_all_LDPruned_mafNA'
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/f50b4/step3.2_ApePhylogeny_f50b4_plot_20210301.R '/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke' 'f50b4_pass_bialleic_all_LDPruned_maf05'
Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/f50b4/step3.2_ApePhylogeny_f50b4_plot_20210301.R '/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke' 'f50b4_pass_bialleic_all_LDPruned_maf10'

```

**START local analyses**

# step2: gds to vcf files

```bash
# Mon Mar  1 02:11:26 2021
run_rscript () {
local WORKDIR=${1}
local GDSFILE=${2}
local OUTPREFIX=${3}
local MAJORREF=${4}
local LOG=${LOGDIR}/step2_Prunedgds2vcf_${OUTPREFIX}.log

Rscript --vanilla /Users/linmeixi/Lab/fin_whale/scripts_analyses/PopStructure/f50b4/step2_Prunedgds2vcf_f50b4_20210223.R ${WORKDIR} ${GDSFILE} ${OUTPREFIX} ${MAJORREF} &>${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Prunedgds2vcf Done"
echo -e "[$(date "+%Y-%m-%d %T")] GIT commit id ${COMMITID}" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Prunedgds2vcf Done" >> ${LOG}
}

TODAY=$(date "+%Y%m%d")
set -o pipefail
HOMEDIR='/Users/linmeixi/Lab/fin_whale/scripts_analyses'
WORKDIR='/Users/linmeixi/google_drive/finwhale/analyses/PopStructure/f50b4/Minke'
LOGDIR=${WORKDIR}/logs
COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

cd ${WORKDIR}
mkdir -p ${LOGDIR}

run_rscript ${WORKDIR} 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA.gds' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA' 'TRUE'
run_rscript ${WORKDIR} 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA.gds' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_mafNA' 'FALSE'

run_rscript ${WORKDIR} 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf05.gds' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf05'
run_rscript ${WORKDIR} 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf10.gds' 'JointCalls_f50b4_filterpass_bialleic_all_LDPruned_maf10'

```

# NOTE:

Final settings used:

minor allele frequency cutoff (maf) = 0.10
pairwise tree was used









