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
# *                                   v.1.0.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
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

#Add necessary programs and files
cp $source/astral.4.10.2.jar .
cp -r $source/lib .

#Copy genetree file
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
else
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
fi

#Modify labels in gene tree
sed -i 's/XX/-/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
sed -i 's/YY/_/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick

#Copy bootrapped gene tree files (if tree=RAxML or tree=FastTree and FastTreeBoot=yes)
if [[ $tree =~ "RAxML" ]]; then
	#Make dir for for bootstraped trees
	mkdir boot
	#Copy RAxML bootstraped trees
	if [[ $update =~ "yes" ]]; then
		cp $path/${alnpathselected}${type}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
			cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_${i}.result boot/
		done
	else
		cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/*bootstrap* boot/
	fi
	#Make a list of bootstraped trees
	ls boot/*bootstrap* > bs-files
elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	#Make dir for for bootstraped trees
	mkdir boot
	#Copy FastTree bootstraped trees
	if [[ $update =~ "yes" ]]; then
		cp $path/${alnpathselected}${type}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
			cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_${i}_*.boot.fast.trees boot/
		done
	else
		cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees boot/
	fi
	#Make a list of bootstraped trees
	ls boot/*.boot.fast.trees > bs-files
fi
#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral
else
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral
fi

#Run ASTRAL
java -jar astral.4.10.2.jar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_updated.tre
if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	#Run Astral bootstrap
	java -jar astral.4.10.2.jar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_updated_allbootstraptrees.tre
	#Extract last tree (main species tree based on non-bootstrapped dataset + bootstrap values mapped on it)
	cat Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_updated_allbootstraptrees.tre | tail -n1 > Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_updated_withbootstrap.tre
fi
#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp Astral*.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral
else
	cp Astral*.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral
fi

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07a
fi
