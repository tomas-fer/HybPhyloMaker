#!/bin/bash
#PBS -l walltime=2d
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=8gb
#PBS -l scratch=1gb
#PBS -N HybPipe7e_concatenated_tree
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                    Script 07e - concatenated species tree                    *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree from selected concatenated genes
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPipe4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#or specify another input genes for concatenation below

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
	module add fasttree-2.1.8
	module add python-3.4.1-intel
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07e
	cd workdir07e
fi

#Add necessary scripts and files
cp $source/AMAS.py .
cp $path/71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .

# Make new dir for results
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
# Copy and modify selected FASTA files
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
do
	cp $path/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
	sed -i 's/(/_/g' ${i}_modif${MISSINGPERCENT}.fas
	sed -i 's/)//g' ${i}_modif${MISSINGPERCENT}.fas
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	sed -i 's/_contigs//g' ${i}_modif${MISSINGPERCENT}.fas
	sed -i 's/.fas//g' ${i}_modif${MISSINGPERCENT}.fas
done

#Prepare concatenated dataset and transform it to phylip format
python3 AMAS.py concat -i *.fas -f fasta -d dna -u fasta -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta
python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
#Copy concatenated file to home
cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
#FastTree for concatenated dataset
if [ ! $LOGNAME == "" ]; then
	fasttreemp -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
else
	fasttree -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
fi
#Copy results to home
cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/species_trees/concatenated

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07e
fi
