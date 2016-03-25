#!/bin/bash
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -N HybPipe6_root_trees
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                          Script 06 - Root gene trees                         *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Modify and root trees using Newick Utilities
# (1) root trees with accession specified in variable OUTGROUP
# (2) remove bootstrap values and branch length information from tree files
# Take trees starting with RAxML_bipartitions* (i.e., best ML tree with BS values) from /concatenated_exon_alignments/selected${CUT}RAxML/
# or trees Assembly*.tre from /concatenated_exon_alignments/selected${CUT}FastTree/
# Run first HybPipe4_missingdataremoval.sh with the same $CUT value 
# and HybPipe5a_RAxML_for_selected.sh or HybPipe5b_FastTree_for_selected.sh with the same $CUT value


if [ ! $LOGNAME == "" ]; then
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
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir05
	cd workdir05
fi

#Make dir for result
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees

#Copy trees from home and make multitree file
if [[ $tree =~ "RAxML" ]]; then
	cp $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bipartitions.* .
	cat *bipartitions.* > trees.newick
else
	cp `find $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree -maxdepth 1 ! -name '*boot*'` .
	cat Assembly*.tre > trees.newick
fi

#Root trees with OUTGROUP using newick utilities
nw_reroot trees.newick $OUTGROUP > trees_rooted.newick
#Remove BS values (necessary for ML-EST via STRAW server?)
nw_topology -Ib trees_rooted.newick > trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Modify names: remove unwanted ends of names, replace '-' by XX and '_' by YY, find every XX after a digit and e and add XXXX, replace XXXX by '-' (to preserve numbers in scientific format)
cat trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick | sed 's/_contigs.fas//g' | sed 's/.fas//g' | sed 's/-/XX/g' | sed 's/_/YY/g' | sed -r 's/([0-9]eXX)/\1XX/g' | sed 's/XXXX/-/g' > tmp && mv tmp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Copy multi-tree file to home
cp trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
cp trees.newick $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
cp trees_rooted.newick $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
#This newick gene tree file can be used by HybPipe7a_astral.sh, HybPipe7b_astrid, HybPipe7c_mrl.sh, HybPipe7d_mpest.sh, HybPipe7e_concatenated.sh

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir05
fi
