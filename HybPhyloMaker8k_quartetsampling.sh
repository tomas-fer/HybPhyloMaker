#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=5:mem=24gb:scratch_local=8gb
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
# *                                   v.1.8.0c                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute quartet sampling using the method of Pease et al. (2018)
#https://github.com/FePhyFoFum/quartetsampling

#Takes concatenated alignment in phylip format from /${tree}/species_trees/concatenated or /${tree}/update/species_trees/concatenated
#Takes species tree (Astral, concatenated FastTree or ExaML)
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

# Copy species tree
if [[ $qstree =~ "Astral" ]]; then
	echo -e "Copying Astral species tree...\n"
	if [[ $update =~ "yes" ]]; then
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral/Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	else
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral/Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	fi
	mv Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre tree.tre
elif [[ $qstree =~ "FastTree" ]]; then
	echo -e "Copying FastTree concatenated species tree...\n"
	if [[ $update =~ "yes" ]]; then
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre .
	else
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre .
	fi
	mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre tree.tre
elif [[ $qstree =~ "ExaML" ]]; then
	echo -e "Copying ExaML concatenated species tree...\n"
	if [[ $update =~ "yes" ]]; then
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML/ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	else
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenatedExaML/ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	fi
	mv ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre tree.tre
else
	echo -e "No species tree selected. Exiting...\n"
	exit 3
fi

# Modify labels in species tree
#replace each second occurrence of ' ' by '_'
sed -i.bak 's/ \([^ ]*\) / \1_/g' tree.tre
#replace all (remaining) spaces by '-'
sed -i.bak2 's/ /-/g' tree.tre

# Compute quartet sampling
echo -e "Calculating quartet sampling...\n"
#clone QS GitHub
git clone https://github.com/FePhyFoFum/quartetsampling.git &> quartetsampling.log
#run main QS script
python3 quartetsampling/pysrc/quartet_sampling.py --tree tree.tre --align concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip --reps 100 --threads 4 --lnlike 2 --results-dir results >> quartetsampling.log


#Process quartet sampling results to produce a tree with colour-labelled nodes
#output trees are modified with sed&grep and plotted in R to resemble trees in the publication
echo -e "Plotting tree with quartet sampling scores...\n"
cd results #enter results folder
#1. tree with 'qc' values (used later for coloring nodes)
#removes everything within '[ ]', i.e. df values after species names
sed -e 's/\[[^][]*\]//g' RESULT.labeled.tre.qc > RESULT.labeled.tre.qc.modif
#removes 'qc='
sed -i 's/qc=//g' RESULT.labeled.tre.qc.modif
#2. tree with all values 'qc/qd/qi' (used for plotting the tree and three scores)
#removes everything within [] except 'score=...' and '['
sed 's/\[\&[^][]*\,//g' RESULT.labeled.tre.figtree > RESULT.labeled.tre.figtree.modif
#removes everything within '[ ]', i.e. df values after species names
sed -i 's/\[[^][]*\]//g' RESULT.labeled.tre.figtree.modif
#removes 'score='
sed -i 's/score=//g' RESULT.labeled.tre.figtree.modif
#removes ']'
sed -i 's/\]//g' RESULT.labeled.tre.figtree.modif
#take only tree, i.e., make newick file
grep tree1 RESULT.labeled.tre.figtree.modif | sed 's/^.*=//' > RESULT.labeled.tre.figtree.modif.nwk
#remove '  ' which usually left at the beginning of the tree
sed -i 's/  //' RESULT.labeled.tre.figtree.modif.nwk
#remove 'QS1'
sed -i 's/QS1//' RESULT.labeled.tre.figtree.modif.nwk

#Copy script from source
if [[ $PBS_O_HOST == *".cz" ]]; then
	cp $source/plotQStree.R .
else
	cp ../$source/plotQStree.R .
fi

#Plot tree with colour-labelled nodes
R --slave -f plotQStree.R $OUTGROUP > R.log 2>&1
rm plotQStree.R
cd ..

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
