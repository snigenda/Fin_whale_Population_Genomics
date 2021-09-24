
# step1: call genotypes on four individuals and joint call with the testsix dataset

Redo the WGSproc processes since the predicted insert size was misspecified earlier in Sept 2020.

## WGSproc1-3

Scripts used: 

```bash
# Wed Jan 13 00:20:34 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Preprocessing

bash WGSproc1_qsub_wrapper_baleen.sh "BalMus01" 1 "Minke"
bash WGSproc1_qsub_wrapper_baleen.sh "BalAcu02" 1 "Minke"
bash WGSproc1_qsub_wrapper_baleen.sh "EubGla01" 1 "Minke"
bash WGSproc1_qsub_wrapper_baleen.sh "MegNov01" 1 "Minke"

```

Progress checking:

| SRA        | Name     | WGSproc1a | WGSproc1b | WGSproc1c | WGSproc2 | WGSproc3 | Backup |
| ---------- | -------- | --------- | --------- | --------- | -------- | -------- | ------ |
| SRR1802584 | BalAcu02 | 6162332          |           |           |          | 6197146         |   ALLDONE     |
| SRR5665644 | BalMus01 | 6160677  |  6161462         |  6173621         |  6183886        |          |  ALLDONE    |
| SRR5665640 | EubGla01 | 6163875          |           |           |          | 6196953;6198840         | ALLDONE       |
| SRR5665639 | MegNov01 | 6188607          | 6196543          |           |          | 6214377         |        |

## Some jobs failed, submit manually

```bash
# Sat Jan 16 22:42:09 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Preprocessing
NAME="MegNov01"
REF="Minke"
qsub -l highp,h_rt=72:00:00,h_data=4G,h_vmem=10G -N WGSproc1_b_${NAME}_${REF} WGSproc1_b_Align_baleen.sh ${NAME} ${NAME}_A 1 ${REF}
```
## WGSproc3 resubmit

```bash
###########################################################
# toolbox
WORKSCRIPT=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Preprocessing/WGSproc3_HaplotypeCaller_baleen.sh
# first input index
# second input sample
cleanWGSproc3 () {
local BADID=${1}
local SAMPLE=${2}
cd /u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/preprocessing/${SAMPLE}/Minke/
for ii in ${BADID}; do
rm 06_${SAMPLE}_Minke_MarkDuplicates_HaplotypeCaller_${ii}.log
rm ${SAMPLE}_MarkDuplicates_${ii}.g.vcf.gz
rm ${SAMPLE}_MarkDuplicates_${ii}.g.vcf.gz.tbi
rm WGSproc3_${SAMPLE}_Minke_MarkDuplicates_${ii}_progress.log
done
}

# check if all the jobs were good
checkWGSproc3 () {
local NAME=${1}
local REF=${2}
local NCOUNT=${3}
local LOGFILE=/u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/baleen_WGSproc3_20210112.out.txt
local PCOUNT=$(grep "Done WGSproc3_HaplotypeCaller for ${NAME} ${REF}" ${LOGFILE} | wc -l)

if [[ ${PCOUNT} -eq ${NCOUNT} ]]; then
echo "All success"
else
echo "Fail"
fi
}

###########################################################
# Sat Jan 16 23:14:14 2021
cleanWGSproc3 "07 08 09 10 11 16 17" "BalMus01"
NAME="BalMus01"
REF="Minke"
qsub -t 7-11 -N WGSproc3_${NAME}_${REF} ${WORKSCRIPT} ${NAME} ${REF}
qsub -t 16-17 -N WGSproc3_${NAME}_${REF} ${WORKSCRIPT} ${NAME} ${REF}

###########################################################
# Sun Jan 17 12:01:52 2021
cleanWGSproc3 "21 30 63 64 78 82 92" "BalMus01"
BADID="21 30 63 64 78 82 92"
NAME="BalMus01"
REF="Minke"
for ii in ${BADID}; do
echo $ii
sleep 5
qsub -t $ii -l h_rt=23:59:00,h_data=20G,h_vmem=24G -N WGSproc3_${NAME}_${REF} ${WORKSCRIPT} ${NAME} ${REF}
done

###########################################################
# Sun Jan 17 12:14:22 2021
# speed up the EubGla01, move to non-highp nodes
NAME="EubGla01"
REF="Minke"
qsub -t 21-96 -l h_rt=23:59:00,h_data=20G,h_vmem=24G -N WGSproc3_${NAME}_${REF} ${WORKSCRIPT} ${NAME} ${REF}

###########################################################
# Mon Jan 18 13:48:50 2021
NAME="BalMus01"
REF="Minke"
cleanWGSproc3 "16" "BalMus01"

qsub -t 16 -N WGSproc3_${NAME}_${REF} ${WORKSCRIPT} ${NAME} ${REF}

###########################################################
# Tue Jan 19 00:35:39 2021
checkWGSproc3 "EubGla01" "Minke" "96"
```

## WGSproc4

```bash
# Tue Jan 19 11:29:25 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Preprocessing

qsub -N WGSproc4_BalAcu02_Minke WGSproc4_FinalCheck_baleen.sh "BalAcu02" "Minke"
qsub -N WGSproc4_BalMus01_Minke WGSproc4_FinalCheck_baleen.sh "BalMus01" "Minke"
qsub -N WGSproc4_EubGla01_Minke WGSproc4_FinalCheck_baleen.sh "EubGla01" "Minke"
qsub -N WGSproc4_MegNov01_Minke WGSproc4_FinalCheck_baleen.sh "MegNov01" "Minke"

```

# TODO: step2: more completed datasets


