#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=24gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8h_neighbour_network
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8h_neighbour_network
#$ -o HybPhyloMaker8h_neighbour_network.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                     Script 08h - neighbour network in R                      *
# *                                   v.1.6.6                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute neighbour network from selected concatenated genes using phangorn package in R
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6b_FastTree_for_selected.sh to create the concatenated dataset of selected genes

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8h is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add R-3.4.3-gcc
	module add debian9-compat
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8h is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08h
	cd workdir08h
	#Add necessary modules
	module load tools/R/3.4.1
else
	echo -e "\nHybPhyloMaker8h is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08h
	cd workdir08h
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
	echo -en "...with corrected reading frame"
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
fi

#Check necessary file
echo -ne "\nTesting if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta" ]; then
		echo -e "OK\n"
	else
		echo -e "no concatenated alignment file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/'. Run first HybPhyloMaker8e_concatenatedFastTree.sh. Exiting..."
		rm -d ../workdir08h 2>/dev/null
		exit 3
	fi
else
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta" ]; then
		echo -e "OK\n"
	else
		echo -e "no concatenated alignment file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/'. Run first HybPhyloMaker8e_concatenatedFastTree.sh. Exiting..."
		rm -d ../workdir08h 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/NeighbourNetwork" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/NeighbourNetwork' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08h 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/NeighbourNetwork" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/NeighbourNetwork' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08h 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08h)" ]; then
		echo -e "Directory 'workdir08h' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08h 2>/dev/null
		exit 3
	fi
fi

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/NeighbourNetwork
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/NeighbourNetwork
fi

# Copy concatenated alignment
echo -e "Copying concatenated alignment...\n"
if [[ $update =~ "yes" ]]; then
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta .
else
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta .
fi


#Compute NeighbourNetwork using phangorn
echo -e "Computing NeighbourNetwork for concatenated dataset...\n"
mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta file.fasta
R -q -e "library(phangorn);a<-read.phyDat('file.fasta',format='fasta',type='DNA');m<-dist.hamming(a);nnet<-neighborNet(m);write.nexus.networx(nnet,file='n.nex')"
mv n.nex NeighbourNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex

#Removing '_cpDNA' from names in network
sed -i.bak 's/_cpDNA//g' NeighbourNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp NeighbourNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/NeighbourNetwork
else
	cp NeighbourNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/NeighbourNetwork
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08h
fi

echo -e "\nHybPhyloMaker8h finished...\n"
