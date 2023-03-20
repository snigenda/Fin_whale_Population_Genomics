#!/bin/bash
#$ -l h_data=10G,h_rt=20:00:00
#$ -wd <homedir>
#$ -o <homedir>/reports/baleen_make_seqdict_20200910.out.txt
#$ -e <homedir>/reports/baleen_make_seqdict_20200910.err.txt
#$ -m abe

# reference: https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk


conda activate gentools

set -eo pipefail

# define variables
HOMEDIR=<homedir2>/finwhale
WORK=${1}
REFERENCE_SEQ=${2} # should be name of reference seq, extension should be fasta but not fna
LOG=${WORK}/${REFERENCE_SEQ/fasta/log}

echo -e "[$(date "+%Y-%m-%d %T")] Start making refseq dictionary in ${WORK} for reference sequence ${REFERENCE_SEQ}"
echo -e "[$(date "+%Y-%m-%d %T")] Start making refseq dictionary in ${WORK} for reference sequence ${REFERENCE_SEQ}" > ${LOG}


# generate index for reference sequence
cd $WORK
bwa index ${REFERENCE_SEQ} &>> ${LOG}

samtools faidx ${REFERENCE_SEQ} &>> ${LOG}

# generate dictionary (step 4c)
picard -Xmx4G CreateSequenceDictionary \
REFERENCE=${REFERENCE_SEQ} \
OUTPUT=${REFERENCE_SEQ/fasta/dict} &>> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] Done"
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
conda deactivate