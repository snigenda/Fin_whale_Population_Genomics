# Generate SIFT database for minke whale genome

# Database creation pipeline: https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB
# SIFT4G: https://github.com/rvaser/sift4g
# Note: Used modified code as detailed in https://github.com/rvaser/sift4g/issues/10

# Used uniref90.fasta downloaded May 1 2020 (109653977 sequences)

### Set variables, paths
PARENT_DIR=~/project/programs/sift/SIFT_databases/GCF_000493695.1_BalAcu1.0
ORG=balaenoptera_acutorostrata
ORG_VERSION=GCF_000493695.1_BalAcu1.0

SIFT4G_PATH=~/project/programs/sift/sift4g_mod/bin/sift4g
PROTEIN_DB=~/project/programs/sift/SIFT_databases/uniref90.fasta

### Set up PARENT_DIR
mkdir ${PARENT_DIR}
cd ${PARENT_DIR}
mkdir chr-src
mkdir dbSNP
mkdir gene-annotation-src

### Get genome fasta and gene annotation
scp jarobins@hoffman2.idre.ucla.edu:/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta .
scp jarobins@hoffman2.idre.ucla.edu:/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.gtf .

gzip GCF_000493695.1_BalAcu1.0_genomic.gtf
mv GCF_000493695.1_BalAcu1.0_genomic.gtf.gz gene-annotation-src/

### Split genome into individual fasta records, then gzip each
# faSplit available from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
faSplit byname GCF_000493695.1_BalAcu1.0_genomic.fasta chr-src/
for i in chr-src/*.fa ; do gzip ${i} ; done

### Generate config file
cat > config <<EOF
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=2
MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial

PARENT_DIR=${PARENT_DIR}
ORG=${ORG}
ORG_VERSION=${ORG_VERSION}

#Running SIFT 4G
SIFT4G_PATH=${SIFT4G_PATH}
PROTEIN_DB=${PROTEIN_DB}

# Sub-directories, don't need to change
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=Log.txt
ZLOGFILE=Log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP

# Doesn't need to change
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord

EOF

### Run make-SIFT-db-all.pl
cd ~/project/programs/sift/scripts_to_build_SIFT_db/
CONFIG=${PARENT_DIR}/config
LOG=${PARENT_DIR}/makeSiftDB.log
/usr/bin/time -v perl ./make-SIFT-db-all.pl -config ${CONFIG} |& tee ${LOG}

### Archive
cd ~/project/programs/sift/SIFT_databases
mv GCF_000493695.1_BalAcu1.0 GCF_000493695.1_BalAcu1.0_makeSiftDB
mv GCF_000493695.1_BalAcu1.0_makeSiftDB/GCF_000493695.1_BalAcu1.0 .
tar -czf GCF_000493695.1_BalAcu1.0_makeSiftDB.tar.gz GCF_000493695.1_BalAcu1.0_makeSiftDB
# rm -r GCF_000493695.1_BalAcu1.0_makeSiftDB/

### Variant annotation example
# NOTE: VCF input file must be uncompressed
# NOTE: SIFT skips any sites without '.' or 'PASS' in FILTER column
# NOTE: Names of files in SIFT database must match chromosome names in VCF file, except 
# that database files must be named e.g. 1.gz rather than chr1.gz, even if 'chr1' is the 
# name used in the VCF file (rename database files and/or chromosomes in VCF files as needed)
SIFT=~/project/programs/sift/SIFT4G_Annotator_v2.4.jar
DB=~/project/programs/sift/SIFT_databases/GCF_000493695.1_BalAcu1.0
VCF=myvariants.vcf
LOG=myvariants_SIFT.log
java -Xmx8g -jar ${SIFT} -c -t -r . -d ${DB} -i ${VCF} |& tee ${LOG}

