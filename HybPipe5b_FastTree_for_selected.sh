#!/bin/bash
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -l mem=4gb
#PBS -N HybPipe5b_FastTree_for_selected
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                   Script 05b - FastTree gene tree building                   *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute gene trees using FastTree for selected genes
# Selection is based on maximum missing data per sample allowed ($MISSINGPERCENT) and minimum species percentage presence per assembly ($SPECIESPRESENCE)
# Edit these values in settings.cfg
# Run first HybPipe4_missingdataremoval.sh with the same settings

#Complete path and set configuration for selected location
if [ ! $LOGNAME == "" ]; then
	echo "Metacentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
		#Add necessary modules
	module add fasttree-2.1.8
	module add raxml-8.2.4
	module add perl-5.10.1
else
	echo "Local..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir05b
	cd workdir05b
fi

#Add necessary scripts and files
cp $source/catfasta2phyml.pl .
cp $source/CompareToBootstrap.pl .
cp $source/MOTree.pm .
cp $path/71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
# Copy and modify selected FASTA files
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
do
	cp $path/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in FastTree)
	sed -i 's/(/_/g' ${i}_modif${MISSINGPERCENT}.fas
	sed -i 's/)//g' ${i}_modif${MISSINGPERCENT}.fas
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	sed -i 's/_contigs//g' ${i}_modif${MISSINGPERCENT}.fas
	sed -i 's/.fas//g' ${i}_modif${MISSINGPERCENT}.fas
done
#Make a list of all fasta files
ls *.fas | cut -d"." -f1 > FileForFastTree.txt
#Make dir for results
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree

for file in $(cat FileForFastTree.txt)
do
	#FastTree
	if [ ! $LOGNAME == "" ]; then
		fasttreemp -nt ${file}.fas > ${file}.fast.tre
	else
		fasttree -nt ${file}.fas > ${file}.fast.tre
	fi
	cp *$file.fast.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
done

#Bootstrap using FastTree 
if [[ $FastTreeBoot =~ "yes" ]]; then
	for file in $(cat FileForFastTree.txt)
	do
		#Generate 100 replicated datasets using RAxML -f j
		raxmlHPC -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS
		#Loop over replicates and calculate FastTree for each of them
		for i in {0..99}
		do
			fasttreemp -nt ${file}.fas.BS${i} > ${file}.BS${i}.fast.tre
		done
		#Combine all bootstrap trees to a single file
		cat ${file}.BS*.fast.tre > ${file}.boot.fast.trees
		#Delete BS files and trees
		rm *.BS*
		rm ${file}.BS*.fast.tre
		#Map bootstrap support values onto the original tree
		perl ./CompareToBootstrap.pl -tree *${file}.fast.tre -boot ${file}.boot.fast.trees > ${file}.boot.fast.tre
		#Copy bootstrap trees and a final tree with bootstrap values to home
		cp ${file}.boot.fast.trees $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
		cp ${file}.boot.fast.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	done
fi

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir05b
fi

