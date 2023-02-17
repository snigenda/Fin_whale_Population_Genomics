# edits from previous versions

1. changed the population sizes to match the supplemental tables.
2. added tracking statistics for important generations.

> notes: in population split, the size change is implemented at the time of the change. but in population size change, the size change is implemented at the end of the change.

# submit the FINAL SLiM scripts

```bash
# 2022-03-16 11:53:32
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/revisions_SLiM_final/

submit_slim() {
    local SCRIPTNAME=${1}
    qsub -N ${SCRIPTNAME} wrapper_revisions_SLiM_final_20220316.sh "${SCRIPTNAME}" "/u/project/rwayne/meixilin/fin_whale/analyses/scripts/revisions_SLiM_final/${SCRIPTNAME}.slim"
}

submit_slim 'finwhale_2pop_AncChange_AsymMig_20220316' # 2025259.1-25:1
submit_slim 'finwhale_2pop_AncChange_NoMig_20220316' # 2025260.1-25:1
submit_slim 'finwhale_ENP_3Epoch_20220316' # 2025262.1-25:1
submit_slim 'finwhale_ENP_3Epoch_recovery_20220316' # 2025263.1-25:1
```

Resubmit slim for some jobs that failed the memory limits

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/revisions_SLiM_final/

# 2022-03-19 15:25:55
resubmit_slim() {
    local SCRIPTNAME=${1}
    local ID=${2}
    qsub -t ${ID} -N ${SCRIPTNAME}.${ID} wrapper_revisions_SLiM_final_20220316.sh "${SCRIPTNAME}" "/u/project/rwayne/meixilin/fin_whale/analyses/scripts/revisions_SLiM_final/${SCRIPTNAME}.slim"
}

resubmit_slim 'finwhale_2pop_AncChange_AsymMig_20220316' '6' # 2059806.6-6:1
```

# original scripts folder

/Users/linmeixi/Lab/fin_whale/scripts_analyses/revisions_SLiM_final

