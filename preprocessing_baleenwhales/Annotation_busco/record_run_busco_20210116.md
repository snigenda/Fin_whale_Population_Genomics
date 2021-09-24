# step1.0: run busco

```bash
# testing runs
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

REF="Minke"
BUSCODB="cetartiodactyla_odb10"
qsub baleen_genomes/Annotation_busco/step1_run_busco_20210116.sh ${REF} ${BUSCODB}
# 6192472
```

# step1: run in larger scales

```bash
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts
REFLIST=("Minke" "Bryde" "Fin" "Blue" "Humpback" "Minke2" "Blue2")
# BUSCODBLIST=("cetartiodactyla_odb10" "laurasiatheria_odb10" "mammalia_odb10" "tetrapoda_odb10")
BUSCODBLIST=("cetartiodactyla_odb10" "mammalia_odb10")
REFLIST=("Minke" "Bryde" "Blue")

for BUSCODB in ${BUSCODBLIST[@]}; do
for REF in ${REFLIST[@]}; do
sleep 5
qsub -N busco_${REF}_${BUSCODB:0:5} baleen_genomes/Annotation_busco/step1_run_busco_20210116.sh ${REF} ${BUSCODB}
done
done
```

## Bryde's whale had some trouble in transcript format

```bash
# Sun Jan 17 10:49:19 2021
cd /u/project/rwayne/snigenda/finwhale/cetacean_genomes/brydes_whale_genome/DNAzoo_data/Balaenoptera_edeni/genome/maker
WORKSCRIPT=/u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Annotation_busco/convert_gap2N.py
FILE="Balaenoptera_edeni_Balaenoptera_edeni_HiC.fasta_v2.functional.transcripts.fasta"
python ${WORKSCRIPT} --infasta ${FILE} --outfasta ${FILE:0:73}2N.fasta
echo $?
```

## Resubmit Bryde's whale

```bash
resub_busco() {
local REF=${1}
cd /u/project/rwayne/meixilin/fin_whale/analyses/reports/baleen_genomes/
rm BUSCO_${REF}*.log
cd /u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/busco
rm -r ${REF}*/
echo $?
}

# Sun Jan 17 12:04:55 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

BUSCODBLIST=("cetartiodactyla_odb10" "mammalia_odb10")
REFLIST=("Bryde" "Fin" "Humpback" "Minke2" "Blue2")
for BUSCODB in ${BUSCODBLIST[@]}; do
for REF in ${REFLIST[@]}; do
sleep 10
qsub -N busco_${REF}_${BUSCODB:0:5} baleen_genomes/Annotation_busco/step1_run_busco_20210116.sh ${REF} ${BUSCODB}
done
done
```

## Bug in the original script about proteinmode

```bash
resub_busco "Minke2"
resub_busco "Blue2"

cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts

REFLIST=("Minke2" "Blue2")
# BUSCODBLIST=("cetartiodactyla_odb10" "laurasiatheria_odb10" "mammalia_odb10" "tetrapoda_odb10")
BUSCODBLIST=("cetartiodactyla_odb10" "mammalia_odb10")

for BUSCODB in ${BUSCODBLIST[@]}; do
for REF in ${REFLIST[@]}; do
sleep 10
qsub -N busco_${REF}_${BUSCODB:0:5} baleen_genomes/Annotation_busco/step1_run_busco_20210116.sh ${REF} ${BUSCODB}
done
done
```

## Plot the busco analyses

1. The "rna.fna" (Minke) and "rna_from_genomic.fna" (Minke2) output is the same.

```bash
# Tue Jan 19 00:11:53 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Annotation_busco

bash step2_plot_busco_20210116.sh

```

## Archive the output

```bash
# Tue Jan 19 00:27:42 2021
cd /u/project/rwayne/meixilin/fin_whale/analyses/scripts/baleen_genomes/Annotation_busco

qsub step3_archive_busco_20210119.sh # JOB ID: 6217241

# Tue Jan 19 10:41:08 2021
cd /u/scratch/m/meixilin/finwhale/analyses/baleen_genomes/busco
WORKDIR=/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/busco/

rsync -ahv busco_output_20210119.md5sum ${WORKDIR}
rsync -ahv *.tar.gz ${WORKDIR}
```

# Final output:

See this folder:
`/u/project/rwayne/meixilin/fin_whale/analyses/baleen_genomes/busco/BUSCO_summaries_20210118`

