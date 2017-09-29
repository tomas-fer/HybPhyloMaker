#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker7_root_trees
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker7_root_trees
#$ -o HybPhyloMaker7_root_trees.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                          Script 07 - Root gene trees                         *
# *                                   v.1.5.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2017 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Modify and root trees using Newick Utilities
# (1) root trees with accession specified in variable OUTGROUP
# (2) remove bootstrap values and branch length information from tree files
# Take trees starting with RAxML_bipartitions* (i.e., best ML tree with BS values) from /concatenated_exon_alignments/selected${CUT}RAxML/
# or trees Assembly*.tre from /concatenated_exon_alignments/selected${CUT}FastTree/
# Run first HybPhyloMaker5_missingdataremoval.sh with the same $CUT value 
# and HybPhyloMaker6a_RAxML_for_selected.sh or HybPhyloMaker6b_FastTree_for_selected.sh with the same $CUT value


if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker7 is running on MetaCentrum..."
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
	module add newick-utils-13042016
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker7 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir07
	cd workdir07
	#Add necessary modules
	module load bioinformatics/newickutilities/0.0
else
	echo -e "\nHybPhyloMaker7 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir07
	cd workdir07
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection\n"
else
	echo -e "\n"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update" ]; then
			if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/" ]; then
				if [ "$(ls -A $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree})" ]; then
					echo -e "OK\n"
				else
					echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}' is empty. Exiting...\n"
					rm -d ../workdir07 2>/dev/null
					exit 3
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}' is missing. Exiting...\n"
				rm -d ../workdir07 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update' is missing. Exiting...\n"
			rm -d ../workdir07 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
		rm -d ../workdir07 2>/dev/null
		exit 3
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}" ]; then
		if [ "$(ls -A $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree})" ]; then
			echo -e "OK\n"
		else
			echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}' is empty. Exiting...\n"
			rm -d ../workdir07 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}' is missing. Exiting...\n"
		rm -d ../workdir07 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/*.newick 2>/dev/null | wc -w) ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees' already includes some *.newick trees. You are probably going to owerwrite previous results. exists. Delete or remove all *.newick files before running this script again. Exiting...\n"
		rm -d ../workdir07 2>/dev/null
		exit 3
	fi
else
	if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/*.newick 2>/dev/null | wc -w) ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees' already includes some *.newick trees. You are probably going to owerwrite previous results. Delete or remove all *.newick files before running this script again. Exiting...\n"
		rm -d ../workdir07 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir07)" ]; then
		echo -e "Directory 'workdir07' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir07 2>/dev/null
		exit 3
	fi
fi
#If working with updated tree list select specific trees (otherwise copy all trees)
if [[ $update =~ "yes" ]]; then
	cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
	#Make dir for result
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
	#Loop over a list of selected trees
	for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
		if [[ $tree =~ "RAxML" ]]; then
			#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
			if [[ $corrected =~ "yes" ]]; then
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/RAxML_bipartitions.CorrectedAssembly_${i}_modif${MISSINGPERCENT}.result .
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/RAxML_bipartitions.Assembly_${i}_modif${MISSINGPERCENT}.result .
			fi
		else
			#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
			if [[ $corrected =~ "yes" ]]; then
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/CorrectedAssembly_${i}_modif${MISSINGPERCENT}*.tre .
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/Assembly_${i}_modif${MISSINGPERCENT}*.tre .
			fi
			#Delete all bootstrap/ no bootstrap trees
			if [[ $FastTreeBoot =~ "yes" ]]; then
				rm *_modif${MISSINGPERCENT}.fast.tre
			else
				rm *boot*
			fi
		fi
	done
	echo -e "Combining trees..."
	cat *_modif* > trees.newick
else
	#Make dir for result
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	#Copy trees from home and make multitree file
	if [[ $tree =~ "RAxML" ]]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bipartitions.* .
		echo -e "Combining trees..."
		cat *bipartitions.* > trees.newick
	else
		if [[ $FastTreeBoot =~ "yes" ]]; then
			cp `find $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree -maxdepth 1 -name '*boot*'` . 2>/dev/null
		else
			cp `find $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree -maxdepth 1 ! -name '*boot*'` . 2>/dev/null
		fi
		echo -e "Combining trees..."
		cat *.tre > trees.newick
	fi
fi

if [[ $cp =~ "yes" ]]; then
	#Removing '_cpDNA' from gene trees in trees.newick
	sed -i.bak 's/_cpDNA//g' trees.newick
fi

#Root trees with OUTGROUP using newick utilities
if [ -z "$OUTGROUP" ]; then
	echo -e "\nTrees will not be rooted, no outgroup was specified..."
else
	echo -e "\nRerooting trees..."
	nw_reroot -s trees.newick $OUTGROUP > trees_rooted.newick 2>root_log.txt
	nrrooted=$(cat root_log.txt | wc -l)
	if [ $nrrooted -eq 0 ]; then
		echo -e "All trees were rooted"
	else
		echo -e "$nrrooted trees were not rooted (Outgroup: ${OUTGROUP} was not found.)\nHowever, the file trees_rooted.newick is still suitable for ASTRAL, ASTRID and MRL analysis.\nIt will not work for any analysis requiring rooted trees, i.e., MP-EST."
	fi
fi

#Remove BS values (necessary for MP-EST via STRAW server?)
echo -e "\nRemoving BS values from trees..."
if [ -z "$OUTGROUP" ]; then
	nw_topology -Ib trees.newick > trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick
else
	nw_topology -Ib trees_rooted.newick > trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
fi
#Modify names: remove unwanted ends of names, replace '-' by XX and '_' by YY, find every XX after a digit and e and add XXXX, replace XXXX by '-' (to preserve numbers in scientific format)
if [ -z "$OUTGROUP" ]; then
	cat trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick | sed 's/_contigs.fas//g' | sed 's/.fas//g' | sed 's/-/XX/g' | sed 's/_/YY/g' | sed -r 's/([0-9]eXX)/\1XX/g' | sed 's/XXXX/-/g' > tmp && mv tmp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick
else
	cat trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick | sed 's/_contigs.fas//g' | sed 's/.fas//g' | sed 's/-/XX/g' | sed 's/_/YY/g' | sed -r 's/([0-9]eXX)/\1XX/g' | sed 's/XXXX/-/g' > tmp && mv tmp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
fi
#Copy multi-tree file to home
if [[ $update =~ "yes" ]]; then
	cp trees*.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
else
	cp trees*.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
fi
#This newick gene tree file can be used by HybPhyloMaker7a_astral.sh, HybPhyloMaker7b_astrid.sh, and HybPhyloMaker7c_mrl.sh

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir07
fi

echo -e "\nScript HybPhyloMaker7 finished...\n"
