#!/bin/bash
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=8gb
#PBS -N HybPipe7b_Astrid
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                       Script 07b - Astrid species tree                       *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using Astrid methods from trees saved in single gene tree file (with *.newick suffix)
#Take trees from 72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Run first
#(1) HybPipe4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPipe5a_RAxML_for_selected.sh or HybPipe5b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPipe6_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#or specify another input trees below


#Complete path and set configuration for selected location
if [ ! $LOGNAME == "" ]; then
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add astrid-1.0
	#module add python-2.7.6-gcc
	#module add python27-modules-gcc
	#module add jdk-7
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir06b
	cd workdir06b
fi

#Copy genetree file
cp $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
#Modify labels in gene tree
sed -i 's/XX/-/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
sed -i 's/YY/_/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick

#Copy bootrapped gene tree files (if tree=RAxML or tree=FastTree and FastTreeBoot=yes)
if [[ $tree =~ "RAxML" ]]; then
	#Make dir for for bootstraped trees
	mkdir boot
	#Copy RAxML bootstraped trees
	cp $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/*bootstrap* boot/
	#Make a list of bootstraped trees
	ls boot/*bootstrap* > bs-files
elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	#Make dir for for bootstraped trees
	mkdir boot
	#Copy FastTree bootstraped trees
	cp $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees boot/
	#Make a list of bootstraped trees
	ls boot/*.boot.fast.trees > bs-files
fi
#Make dir for results
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid

#Run ASTRID
ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
#ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -m bionj -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	#Run Astral bootstrap
	ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre
	#Extract last tree (main species tree based on non-bootstrapped dataset + bootstrap values mapped on it)
	cat Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre | tail -n1 > Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
fi

#Copy results to home
cp Astrid*.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07b
fi
