

# Gather raw data from hoffman2

```bash
cd /u/project/rwayne/pnunez/Results/Demography/fsc2/1DModels/resultsSummaries

LOCAL=~/Lab/finwhale_manuscript/data/demography/fsc26/estimate/raw_data/
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/pnunez/Results/Demography/fsc2/1DModels/resultsSummaries
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/pnunez/Results/Demography/fsc2/2DModels/resultsSummaries


rsync -ahv --update -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}
```