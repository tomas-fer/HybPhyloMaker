#!/bin/bash
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=12
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=1gb
#PBS -N MRL
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                        Script 07c - MRL species tree                         *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using MRL methods using RAxML from trees saved in single gene tree file (with *.newick suffix)
#Take trees from /concatenated_exon_alignments/selected${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Run first
#(1) HybPipe4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPipe5a_RAxML_for_selected.sh or HybPipe5b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPipe6_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#or specify another input trees below

#Complete path and set configuration for selected location
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
	module add jdk-1.6.0
	module add raxml-8.2.4
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07c
	cd workdir07c
fi

#Add necessary programs and files
cp $source/mrp.jar .

#Copy genetree file
cp $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .

#Make dir for results
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/MRL

#Make MRP matrix
#java -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRPmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex NEXUS
#Make MRL matrix
java -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip PHYLIP -randomize

#Make 100 fast bootstrap ML trees using RAxML
raxmlHPC -T 12 -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -# 100 
#Modify labels in RAxML bipartitions (XX and YY to ' ')
sed -i 's/XX/ /g' RAxML_bipartitions.MRLresult
sed -i 's/YY/ /g' RAxML_bipartitions.MRLresult

#Copy results to home
cp *MR* $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/MRL

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07c
fi
