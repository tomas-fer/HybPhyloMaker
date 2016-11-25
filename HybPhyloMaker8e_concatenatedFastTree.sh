#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2d
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -l mem=24gb
#PBS -l scratch=8gb
#PBS -N HybPhyloMaker7e_concatenated_tree
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker7e_concatenated_tree
#$ -o HybPhyloMaker7e_concatenated_tree.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                    Script 07e - concatenated species tree                    *
# *                                   v.1.3.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree from selected concatenated genes
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPhyloMaker4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker7e is running on MetaCentrum..."
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
	module add newick-utils-1.6
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker7e is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir07e
	cd workdir07e
	#Add necessary modules
	module load bioinformatics/fasttree/2.1.8
	module load bioinformatics/anaconda3/2.3.0
	module load bioinformatics/newickutilities/0.0
else
	echo -e "\nHybPhyloMaker7e is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir07e
	cd workdir07e
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection\n"
else
	echo -e "\n"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
		if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
			echo -e "OK\n"
		else
			echo -e "no alignmenet files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
			rm -d ../workdir07e 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
		rm -d ../workdir07e 2>/dev/null
		exit 3
	fi
else
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
		if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
			echo -e "OK\n"
		else
			echo -e "no alignmenet files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
			rm -d ../workdir07e 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
		rm -d ../workdir07e 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir07e 2>/dev/null
		exit 3
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir07e 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir07e)" ]; then
		echo -e "Directory 'workdir07e' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir07e 2>/dev/null
		exit 3
	fi
fi

#Add necessary scripts and files
cp $source/AMAS.py .
#Copy list of genes
if [[ $update =~ "yes" ]]; then
	cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
	mv selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
else
	cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
fi

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

# Copy and modify selected FASTA files
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d"_" -f2); do
	#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
	cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_${i}_modif${MISSINGPERCENT}.fas .
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
	sed -i.bak 's/(/_/g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	sed -i.bak 's/)//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	sed -i.bak 's/_contigs//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	sed -i.bak 's/.fas//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
done

#Prepare concatenated dataset and transform it to phylip format
echo -e "Concatenating...\n"
python3 AMAS.py concat -i *.fas -f fasta -d dna -u fasta -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta >/dev/null
python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip >/dev/null
#Copy concatenated file to home
if [[ $update =~ "yes" ]]; then
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

#FastTree for concatenated dataset
echo -e "Computing FastTree for concatenated dataset...\n"
if [[ $location == "1" ]]; then
	fasttreemp -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
elif [[ $location == "2" ]]; then
	FastTreeMP -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
else
	$fasttreebin -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
fi

#Removing '_cpDNA' from names in alignment
sed -i.bak 's/_cpDNA//g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre

#(Re)root a final concatenated species tree with $OUTGROUP
if [ -n "$OUTGROUP" ]; then
	nw_reroot concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $OUTGROUP > tmp && mv tmp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
fi

#Modify labels in concatenated trees
sed -i.bak2 's/-/ /g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
sed -i.bak2 's/_/ /g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir07e
fi

echo -e "\nHybPhyloMaker7e finished...\n"
