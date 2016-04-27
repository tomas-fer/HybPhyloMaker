#!/bin/bash
#PBS -l walltime=2d
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -l mem=24gb
#PBS -l scratch=8gb
#PBS -N HybPipe7e_concatenated_tree
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                    Script 07e - concatenated species tree                    *
# *                                   v.1.0.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
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

#Add necessary scripts and files
cp $source/AMAS.py .
#Copy list of genes
if [[ $update =~ "yes" ]]; then
	cp $path/${alnpathselected}${type}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
	mv selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
else
	cp $path/${alnpathselected}${type}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
fi

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

# Copy and modify selected FASTA files
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d"_" -f2); do
	#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
	cp $path/${alnpathselected}${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_${i}_modif${MISSINGPERCENT}.fas .
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
	sed -i 's/(/_/g' *ssembly_${i}_modif${MISSINGPERCENT}.fas .
	sed -i 's/)//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas .
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	sed -i 's/_contigs//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas .
	sed -i 's/.fas//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas .
done

#Prepare concatenated dataset and transform it to phylip format
python3 AMAS.py concat -i *.fas -f fasta -d dna -u fasta -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta
python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
#Copy concatenated file to home
#Make new dir for results
if [[ $update =~ "yes" ]]; then
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

#FastTree for concatenated dataset
if [ ! $LOGNAME == "" ]; then
	fasttreemp -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
else
	fasttree -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
fi
#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/update/species_trees/concatenated
else
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/species_trees/concatenated
fi

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07e
fi
