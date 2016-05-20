#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -N HybPipe6_root_trees
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe6_root_trees
#$ -o HybPipe6_root_trees.log

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                          Script 06 - Root gene trees                         *
# *                                   v.1.0.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Modify and root trees using Newick Utilities
# (1) root trees with accession specified in variable OUTGROUP
# (2) remove bootstrap values and branch length information from tree files
# Take trees starting with RAxML_bipartitions* (i.e., best ML tree with BS values) from /concatenated_exon_alignments/selected${CUT}RAxML/
# or trees Assembly*.tre from /concatenated_exon_alignments/selected${CUT}FastTree/
# Run first HybPipe4_missingdataremoval.sh with the same $CUT value 
# and HybPipe5a_RAxML_for_selected.sh or HybPipe5b_FastTree_for_selected.sh with the same $CUT value


if [[ $PBS_O_HOST == *".cz" ]]; then
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
	module add newick-utils-1.6
elif [[ $HOSTNAME == *local* ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir06
	cd workdir06
	#Add necessary modules
	module load bioinformatics/newickutilities/0.0
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir06
	cd workdir06
fi
#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=80concatenated_exon_alignments_corrected
	alnpathselected=81selected_corrected
	treepath=82trees_corrected
else
	alnpath=70concatenated_exon_alignments
	alnpathselected=71selected
	treepath=72trees
fi
#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	type="_cp"
else
	type=""
fi

#If working with updated tree list select specific trees (otherwise copy all trees)
if [[ $update =~ "yes" ]]; then
	cp $path/${alnpathselected}${type}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
	#Make dir for result
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
	#Loop over a list of selected trees
	for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
		#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
		if [[ $corrected =~ "yes" ]]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/CorrectedAssembly_${i}_modif*.tre .
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/Assembly_${i}_modif*.tre .
		fi
	done
	cat *.tre > trees.newick
else
	#Make dir for result
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	#Copy trees from home and make multitree file
	if [[ $tree =~ "RAxML" ]]; then
		cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bipartitions.* .
		cat *bipartitions.* > trees.newick
	else
		cp `find $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree -maxdepth 1 ! -name '*boot*'` .
		cat *.tre > trees.newick
	fi
fi

#Root trees with OUTGROUP using newick utilities
nw_reroot trees.newick $OUTGROUP > trees_rooted.newick
#Remove BS values (necessary for ML-EST via STRAW server?)
nw_topology -Ib trees_rooted.newick > trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Modify names: remove unwanted ends of names, replace '-' by XX and '_' by YY, find every XX after a digit and e and add XXXX, replace XXXX by '-' (to preserve numbers in scientific format)
cat trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick | sed 's/_contigs.fas//g' | sed 's/.fas//g' | sed 's/-/XX/g' | sed 's/_/YY/g' | sed -r 's/([0-9]eXX)/\1XX/g' | sed 's/XXXX/-/g' > tmp && mv tmp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Copy multi-tree file to home
if [[ $update =~ "yes" ]]; then
	cp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
	cp trees.newick $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
	cp trees_rooted.newick $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
else
	cp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	cp trees.newick $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	cp trees_rooted.newick $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
fi
#This newick gene tree file can be used by HybPipe7a_astral.sh, HybPipe7b_astrid.sh, HybPipe7c_mrl.sh, HybPipe7d_mpest.sh, HybPipe7e_concatenated.sh

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir06
fi
