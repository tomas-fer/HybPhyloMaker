#!/bin/bash
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=16gb
#PBS -N HybPipe7a_Astral
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                       Script 07a - Astral species tree                       *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using Astral methods from trees saved in single gene tree file (with *.newick suffix)
#Take trees from 72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
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
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07a
	cd workdir07a
fi

#Add necessary programs and files
cp $source/astral.4.10.2.jar .
cp -r $source/lib .

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
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral

#Run ASTRAL
java -jar astral.4.10.2.jar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	#Run Astral bootstrap
	java -jar astral.4.10.2.jar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre
	#Extract last tree (main species tree based on non-bootstrapped dataset + bootstrap values mapped on it)
	cat Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre | tail -n1 > Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
fi
#Copy results to home
cp Astral*.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07a
fi
