# step1: extract the vcf file

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

for TYPE in 'HIGH' 'MODERATE' 'LOW' 'MODIFIER'; do
echo $TYPE
sleep 2
qsub -N snpEff_impactCDS_${TYPE} get_ALLregions_CDS/snpEff_impactCDS/all50/step1_snpEff_impactCDS_extract_ALLregions_vcf_20210415.sh ${TYPE}
done
```

| Type | Job ID |      |
| ---- | ------ | ---- |
|    HIGH  |    7327758    |      |
|   MODERATE   |  7327759      |      |
|    LOW  |    7327760    |      |
|    MODIFIER  |   7327761     |      |


# transfer files

```bash
# Thu Apr 15 16:52:41 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/get_ALLregions_CDS/all50/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/JointCalls_all50_filterpass_*_ALLregions_all_snpEff.vcf.gz*

rsync -ahv --update --exclude 'JointCalls_all50_filterpass_MODIFIER_ALLregions_all_snpEff.vcf.gz' -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}

# Fri Apr 16 00:17:44 2021
LOCAL=/Users/linmeixi/google_drive/finwhale/analyses/get_ALLregions_CDS/all50/Minke
REMOTE=meixilin@hoffman2.idre.ucla.edu:/u/project/rwayne/meixilin/fin_whale/analyses/get_ALLregions_CDS/all50/Minke/

rsync -ahv --update -e "ssh -i ~/.ssh/id_hoff_rsa" ${REMOTE} ${LOCAL}
```

# note:

Some LOF does not fall within the canonical coding regions but are STILL very important. In this version, some of the high impact sites were not included.

