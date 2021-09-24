# Installing programs & scripts to generate a SIFT database for a new genome

### Make a directory for SIFT
mkdir ~/project/programs/sift
cd ~/project/programs/sift

### Make a directory for databases
mkdir ~/project/programs/sift/databases
cd ~/project/programs/sift/databases

### Obtain and unzip uniref90.fasta (47 Gb uncompressed)
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gzip -d uniref90.fasta.gz

### Obtain and unzip sift4g and SIFT4G_Create_Genomic_DB:
# https://github.com/rvaser/sift4g/archive/master.zip
# https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB/archive/master.zip

### Install prerequisites (other dependencies may be needed depending on system)
sudo apt-get install libswitch-perl
sudo cpan install Bundle::DBI
sudo cpan install Bio::DB::Fasta
sudo cpan install LWP::Simple

### Install sift4g:
# - Modify code in sift4g/vendor/swsharp/swsharp/src/cpu_module.c as described in 
# https://github.com/rvaser/sift4g/issues/10
# - Run 'make cpu' in ~/project/programs/sift/sift4g/vendor/swsharp/swsharp/src
# - Run 'make' in ~/project/programs/sift/sift4g
# - Test that things are working with '~/project/programs/sift/sift4g/bin/sift4g -h'

### Test SIFT4G_Create_Genomic_DB is operational
# Use test examples at https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB


