#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -N Rpackages_setup
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                              Setup R packages                                *
# *                                   v.1.0.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Download working versions of specific packages and install them to personal library
#ape 3.5 (IMPORTANT, v 3.4 is not compatible with HybPipe)
#data.table 1.9.6
#phytools 0.5-20
#seqinr 3.1-3


#Move to scratch
cd $SCRATCHDIR
#Copy file with settings from home and set variables from settings.cfg
cp $PBS_O_WORKDIR/settings.cfg .
. settings.cfg
. /packages/run/modules-2.0/init/bash
#Add necessary modules
module add R-3.2.3-intel

#Get ape 3.5 and other R packages
wget https://cran.r-project.org/src/contrib/ape_3.5.tar.gz 
wget https://cran.r-project.org/src/contrib/data.table_1.9.6.tar.gz
wget https://cran.r-project.org/src/contrib/phytools_0.5-20.tar.gz
wget https://cran.r-project.org/src/contrib/seqinr_3.1-3.tar.gz

#Install packages to personal (writable) library
mkdir -p /storage/$server/home/$LOGNAME/Rpackages
for package in $(ls *.tar.gz); do
	R CMD INSTALL $package -l /storage/$server/home/$LOGNAME/Rpackages
done

#Use this command before running R in HybPipe scripts
export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"

#Delete scratch
rm -rf $SCRATCHDIR/*