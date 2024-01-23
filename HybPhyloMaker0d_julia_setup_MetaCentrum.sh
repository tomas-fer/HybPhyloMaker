#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=1gb
#PBS -j oe
#PBS -N julia_packages_setup
#PBS -m abe

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *               Script 0d - Setup julia packages on Metacentrum                *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2024 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nScript HybPhyloMaker0d is running on Metacentrum..."
else
	echo -e "\nYou are not on Metacentrum. Install the julia packages by hand. Exiting...\n"
	exit 3
fi

#Move to scratch
cd $SCRATCHDIR
#Copy file with settings from home and set variables from settings.cfg
cp $PBS_O_WORKDIR/settings.cfg .
. settings.cfg

#Add necessary modules
module add r/4.0.0-gcc
module add julia

#Use this command before running julia in HybPhyloMaker scripts
export JULIA_DEPOT_PATH=/storage/${server}/home/${LOGNAME}/.julia

cat << EOF | julia
#check installed packages
import Pkg; Pkg.status()
#install packages
Pkg.add("PhyloNetworks")
Pkg.add("PhyloPlots")
Pkg.add("RCall")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Distributed")
Pkg.build("RCall")
using PhyloPlots
Pkg.status()
EOF

#Delete scratch
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
fi

echo -e "\nScript HybPhyloMaker0d finished...\n"
