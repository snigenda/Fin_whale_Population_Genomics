#!/bin/bash
#
# @version 		v0
# @script		.sh
# @description	record the install process
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jan 15 13:23:52 2021
# Use the python3 version
# https://busco.ezlab.org/busco_userguide.html

###########################################################
## import packages
eval "$(/u/project/rwayne/meixilin/miniconda3/bin/conda shell.bash hook)"
conda activate busco

###########################################################
## def functions

###########################################################
## def variables

###########################################################
## main

conda create --name busco

conda activate busco

conda install -c bioconda -c conda-forge busco=4.1.4

busco --list-datasets

                 # - mammalia_odb10
                 #     - eutheria_odb10
                 #         - euarchontoglires_odb10
                 #             - glires_odb10
                 #             - primates_odb10
                 #         - laurasiatheria_odb10
                 #             - carnivora_odb10
                 #             - cetartiodactyla_odb10

# config files

# /u/project/rwayne/meixilin/miniconda3/envs/busco/share/busco/config.ini