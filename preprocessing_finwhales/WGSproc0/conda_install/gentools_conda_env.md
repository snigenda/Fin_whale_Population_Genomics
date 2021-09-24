This is a record on setting up `gentools` conda environment in `/u/project/rwayne/software/finwhale` 

Date: 11/22/2019

# Install miniconda2

```bash
## 1. download the latest miniconda installer 
WORK=/u/project/rwayne/software/finwhale
mkdir -p $WORK
cd $WORK

wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh

## 2. verify the installer hashes 
HASH="383fe7b6c2574e425eee3c65533a5101e68a2d525e66356844a80aa02a556695"
if [ "$(sha256sum Miniconda2-latest-Linux-x86_64.sh | cut -d " " -f 1)" == "$HASH" ]; then 
	echo "Hash verified"
fi
# 383fe7b6c2574e425eee3c65533a5101e68a2d525e66356844a80aa02a556695  Miniconda2-latest-Linux-x86_64.sh

## 3. run installation script 
chmod 755 Miniconda2-latest-Linux-x86_64.sh 
bash Miniconda2-latest-Linux-x86_64.sh 
## check installation logs 
# /u/project/rwayne/software/finwhale/miniconda_installation.log

# Python 2.7.15 | packaged by conda-forge | (default, Feb 28 2019, 04:00:11) 
# [GCC 7.3.0] on linux2
```


# Creating gentools environment and installing software

```bash
## 1. Create environment and do installation
WORK=/u/project/rwayne/software/finwhale
cd $WORK
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda env create --file gentools_env.yml
## check environment setup logs 
# /u/project/rwayne/software/finwhale/conda_gentools_create.log


## 2. Comment out lines 416-432 of lib/R/bin/javareconf in the environment's folder (eg. ~/anaconda2/envs/gentools/lib/R/bin) like so:

cd /u/project/rwayne/software/finwhale/miniconda2/envs/gentools/lib/R/bin
vi javareconf

# if test -n "${R_ARCH}"; then
#   files="${R_HOME}/etc/Makeconf ${R_HOME}/etc/ldpaths ${R_HOME}/etc${R_ARCH}/Makeconf ${R_HOME}/etc${R_ARCH}/ldpaths"
# else
#   files="${R_HOME}/etc/Makeconf ${R_HOME}/etc/ldpaths"
# fi
# for file in $files; do
#     ${SED-sed} -e "s|JAVA =.\{0,\}|JAVA = $JAVA|" -e "s|JAVA_HOME =.\{0,\}|JAVA_HOME = ${JAVA_HOME}|" -e "s|: \${JAVA_HOME=.\{1,\}|: \${JAVA_HOME=${JAVA_HOME}}|" -e "s|: \${R_JAVA_LD_LIBRARY_PATH=.\{1,\}|: \${R_JAVA_LD_LIBRARY_PATH=${JAVA_LD_LIBRARY_PATH_SH}}|" -e "s|JAVA_LIBS =.\{0,\}|JAVA_LIBS = ${JAVA_LIBS}|g" -e "s|JAVA_LD_LIBRARY_PATH =.\{0,\}|JAVA_LD_LIBRARY_PATH = ${JAVA_LD_LIBRARY_PATH}|" -e "s|JAVAC =.\{0,\}|JAVAC = $JAVAC|" -e "s|JAVAH =.\{0,\}|JAVAH = $JAVAH|" -e "s|JAR =.\{0,\}|JAR = $JAR|" -e "s|JAVA_CPPFLAGS =.\{0,\}|JAVA_CPPFLAGS = ${JAVA_CPPFLAGS}|g"  "${file}" > "${file}.new"
#     if test -f "${file}.new"; then
# 	if test -f "${file}"; then
# 	    rm -rf "${file}" >/dev/null 2>&1
# 	fi
# 	mv -f "${file}.new" "${file}"
#     else
# 	echo "*** cannot create ${file}.new~*** Please run as root if required.~" | ${SED-sed} -e 'y/~/\n/' >&2
# 	exit 1
#     fi
# done


## 3. Go to opt/gatk-3.8 in the environment's folder (eg. ~/anaconda2/envs/gentools)
cd /u/project/rwayne/software/finwhale/miniconda2/envs/gentools/opt/gatk-3.8
# Download/transfer/copy GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 to the folder
wget https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef # downloaded online and transfered 
# Run gatk3-register.sh
bash gatk3-register.sh GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
```



# to call gentools

```bash
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools 
```

Date: 02/10/2020
