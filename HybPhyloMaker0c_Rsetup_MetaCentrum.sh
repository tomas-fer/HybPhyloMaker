#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N Rpackages_setup
#PBS -m abe

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                 Script 0c - Setup R packages on Metacentrum                  *
# *                                   v.1.6.6                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2018 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Download working versions of specific packages and install them to personal library
#ape 3.5 (IMPORTANT, v 3.4 is not compatible with HybPhyloMaker)
#data.table 1.9.6
#phytools 0.5-20
#seqinr 3.1-3
#openxlsx 4.0.0

echo -e "\nScript HybPhyloMaker0c is running on Metacentrum..."

#Move to scratch
cd $SCRATCHDIR
#Copy file with settings from home and set variables from settings.cfg
cp $PBS_O_WORKDIR/settings.cfg .
. settings.cfg
. /packages/run/modules-2.0/init/bash
#Add necessary modules
module add R-3.4.3-gcc

#Get ape 3.5 and other R packages
echo -e "\nDownloading R packages..."
wget https://cran.r-project.org/src/contrib/Archive/ape/ape_5.0.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.9.6.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/phytools/phytools_0.5-20.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/seqinr/seqinr_3.1-3.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/openxlsx/openxlsx_4.0.0.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/phangorn/phangorn_2.5.5.tar.gz

#Install packages to personal (writable) library
mkdir -p /storage/$server/home/$LOGNAME/Rpackages
for package in $(ls *.tar.gz); do
	echo -e "\nInstalling $package..."
	R CMD INSTALL $package -l /storage/$server/home/$LOGNAME/Rpackages
done

#Use this command before running R in HybPhyloMaker scripts
export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"

#Delete scratch
rm -rf $SCRATCHDIR/*

echo -e "\nScript HybPhyloMaker0c finished...\n"
