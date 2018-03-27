#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=4:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker10_requisite_collapse
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 2
#$ -q sThM.q
#$ -l mres=8G,h_data=8G,h_vmem=8G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker10_requisite_collapse
#$ -o HybPhyloMaker10_requisite_collapse.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *  Script 10 - Select trees with requisite taxa, collapse unsupported branches *
# *                                   v.1.6.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2018 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using ASTRAL methods from trees saved in single gene tree file (with *.newick suffix)
#Take trees from 72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted.newick
#First, branch in trees that have lower support than required are collapsed
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6a_RAxML_for_selected.sh or HybPhyloMaker6b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPhyloMaker7_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker10_requisite_collapse is running on MetaCentrum..."
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
	module add jdk-1.6.0
	module add newick-utils-1.6
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker10_requisite_collapse is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir10
	cd workdir10
	#Add necessary modules
	module load java/1.7
	module load bioinformatics/anaconda3/2.3.0
	module load bioinformatics/newickutilities/0.0
	module load bioinformatics/p4/ #???
else
	echo -e "\nHybPhyloMaker10_requisite_collapse is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir10
	cd workdir10
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

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection"
else
	if [[ $requisite =~ "no" ]]; then
		echo -e ""
	fi
fi

if [[ $requisite =~ "yes" ]]; then
	echo -e "...and only with trees with requisite taxa present\n"
else
	echo -e "\n"
fi

#Add necessary programs and files
cp $source/TreeCollapseCL4.jar .

if [[ $requisite =~ "no" ]] && [[ $collapse -eq "0" ]]; then
	echo "Nothing to do..."
fi

#Copy genetree file
if [[ $requisite =~ "yes" ]]; then
	if [[ $update =~ "yes" ]]; then
		if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/trees_rooted_with_requisite.newick" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/trees_rooted_with_requisite.newick .
		else
			if [ -z "$OUTGROUP" ]; then
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
				mv trees.newick trees_rooted.newick
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees_rooted.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
			fi
		fi
	else
		if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/trees_rooted_with_requisite.newick" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/trees_rooted_with_requisite.newick .
		else
			if [ -z "$OUTGROUP" ]; then
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
				mv trees.newick trees_rooted.newick
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees_rooted.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
			fi
		fi
	fi
else
	if [[ $update =~ "yes" ]]; then
		if [ -z "$OUTGROUP" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees.newick .
			mv trees.newick trees_rooted.newick
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees_rooted.newick .
		fi
	else
		if [ -z "$OUTGROUP" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees.newick .
			mv trees.newick trees_rooted.newick
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees_rooted.newick .
		fi
	fi
fi

#Select only trees with requisite taxa and copy the file to home
if [[ $requisite =~ "yes" ]]; then
	if [ ! -f "trees_rooted_with_requisite.newick" ]; then
		echo -e "Selecting trees with requisite taxa..."
		grep -E "$requisitetaxa" trees_rooted.newick > trees_rooted_with_requisite.newick
		#Copy selected trees to home
		if [[ $update =~ "yes" ]]; then
			if [ -z "$OUTGROUP" ]; then
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/trees_with_requisite.newick
			else
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
			fi
		else
			if [ -z "$OUTGROUP" ]; then
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/trees_with_requisite.newick
			else
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
			fi
		fi
	fi
fi

#Collapse trees
if [[ ! $collapse -eq "0" ]]; then
	echo -e "Collapsing tree branches with BS below ${collapse}..."
	unset i
	mkdir trees
	#Put all trees one-by-one to folder 'tree'
	if [[ $requisite =~ "yes" ]]; then
		for x in $(<trees_rooted_with_requisite.newick); do
			echo "$x" > trees/tree$((++i)).tre
		done
		java -jar TreeCollapseCL4.jar -b ${collapse} -d trees/ 1>/dev/null
		cat trees/*coll* > trees_with_requisite_collapsed${collapse}.newick
		if [[ $update =~ "yes" ]]; then
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/collapsed${collapse}
			cp trees_with_requisite_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/collapsed${collapse}
		else
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/collapsed${collapse}
			cp trees_with_requisite_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/collapsed${collapse}
		fi
	else
		for x in $(<trees_rooted.newick); do
			echo "$x" > trees/tree$((++i)).tre
		done
		java -jar TreeCollapseCL4.jar -b ${collapse} -d trees/ 1>/dev/null
		cat trees/*coll* > trees_collapsed${collapse}.newick
		if [[ $update =~ "yes" ]]; then
			mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/collapsed${collapse}
			cp trees_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/collapsed${collapse}
		else
			mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/collapsed${collapse}
			cp trees_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/collapsed${collapse}
		fi
	fi
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir10
fi

echo -e "\nScript HybPhyloMaker10_requisite_collapse finished...\n"
