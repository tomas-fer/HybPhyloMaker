#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=24gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8k_quartet_sampling
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8k_quartet_sampling
#$ -o HybPhyloMaker8k_quartet_sampling.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                        Script 08k - quartet sampling                         *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute quartet sampling using the method of Pease et al. (2018)
#https://github.com/FePhyFoFum/quartetsampling

#Takes concatenated alignment in phylip format from /${tree}/species_trees/concatenated or /${tree}/update/species_trees/concatenated
#Takes Astral species tree
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6b_FastTree_for_selected.sh to create the concatenated dataset of selected genes
#(3) HybPhyloMaker7_roottrees.sh and HybPhyloMaker8a_astral.sh to create Astral input tree

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8k is running on MetaCentrum..."
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
	module add raxml-ng-8
	module add debian9-compat
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8k is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08k
	cd workdir08k
	#Add necessary modules
	module load tools/R/3.4.1
else
	echo -e "\nHybPhyloMaker8k is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08k
	cd workdir08k
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

#Settings for selection and (un)corrected reading frame
if [ -z "$selection" ]; then
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/61mafft_corrected
		alnpath=$type/80concatenated_exon_alignments_corrected
		alnpathselected=$type/81selected_corrected
		treepath=$type/82trees_corrected
		echo -en "...with corrected reading frame"
	else
		mafftpath=$type/60mafft
		alnpath=$type/70concatenated_exon_alignments
		alnpathselected=$type/71selected
		treepath=$type/72trees
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/$selection/61mafft_corrected
		alnpath=$type/$selection/80concatenated_exon_alignments_corrected
		alnpathselected=$type/$selection/81selected_corrected
		treepath=$type/$selection/82trees_corrected
		echo -en "...with corrected reading frame...and for selection: $selection"
	else
		mafftpath=$type/$selection/60mafft
		alnpath=$type/$selection/70concatenated_exon_alignments
		alnpathselected=$type/$selection/71selected
		treepath=$type/$selection/72trees
		echo -en "...and for selection: $selection"
	fi
fi

#Check necessary file
echo -ne "\nTesting if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta" ]; then
		echo -e "OK\n"
	else
		echo -e "no concatenated alignment file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/'. Run first HybPhyloMaker8e_concatenatedFastTree.sh. Exiting..."
		rm -d ../workdir08k 2>/dev/null
		exit 3
	fi
else
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta" ]; then
		echo -e "OK\n"
	else
		echo -e "no concatenated alignment file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/'. Run first HybPhyloMaker8e_concatenatedFastTree.sh. Exiting..."
		rm -d ../workdir08k 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/quartetsampling" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/quartetsampling' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08k 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/quartetsampling" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/quartetsampling' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08k 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08k)" ]; then
		echo -e "Directory 'workdir08k' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08k 2>/dev/null
		exit 3
	fi
fi

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/quartetsampling
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/quartetsampling
fi

# Copy concatenated alignment
echo -e "Copying concatenated alignment...\n"
if [[ $update =~ "yes" ]]; then
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
else
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
fi

# Removing '_cpDNA' from names in alignment
sed -i.bak 's/_cpDNA//g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip

# Copy Astral species tree
echo -e "Copying Astral species tree...\n"
if [[ $update =~ "yes" ]]; then
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral/Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
else
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral/Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
fi
# Modify labels in Astral tree
#replace each second occurrence of ' ' by '_'
sed -i.bak 's/ \([^ ]*\) / \1_/g' Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
#replace all (remaining) spaces by '-'
sed -i.bak2 's/ /-/g' Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre

# Compute quartet sampling
echo -e "Calculating quartet sampling...\n"
#clone QS GitHub
git clone https://github.com/FePhyFoFum/quartetsampling.git &> quartetsampling.log
#run main QS script
python3 quartetsampling/pysrc/quartet_sampling.py --tree Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre --align concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip --reps 100 --threads 4 --lnlike 2 --results-dir results >> quartetsampling.log

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp ./results/* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/quartetsampling
	cp quartetsampling.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/quartetsampling
else
	cp ./results/* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/quartetsampling
	cp quartetsampling.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/quartetsampling
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08k
fi

echo -e "\nHybPhyloMaker8k finished...\n"
