#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=8gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8g_BUCKy
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q mThC.q
#$ -l mres=6G,h_data=6G,h_vmem=6G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8g_BUCKy
#$ -o HybPhyloMaker8g_BUCKy.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                    Script 08g - BUCKY concordant analysis                    *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree and CF from selected concatenated genes
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8g is running on MetaCentrum..."
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
	module add bucky-1.4.4
	module add R-3.4.3-gcc
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8g is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08g
	cd workdir08g
	#Add necessary modules
	module load bioinformatics/bucky/1.4.4
	module load tools/R/3.4.1
else
	echo -e "\nHybPhyloMaker8g is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08g
	cd workdir08g
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

#Settings for selection and (un)corrected reading frame
if [ -z "$selection" ]; then
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/61mafft_corrected
		alnpath=$type/80concatenated_exon_alignments_corrected
		alnpathselected=$type/81selected_corrected
		treepath=$type/82trees_corrected
		echo -en "...with corrected reading frame"
	else
		mafftpath=$type/60mafft
		alnpath=$type/70concatenated_exon_alignments
		alnpathselected=$type/71selected
		treepath=$type/72trees
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/$selection/61mafft_corrected
		alnpath=$type/$selection/80concatenated_exon_alignments_corrected
		alnpathselected=$type/$selection/81selected_corrected
		treepath=$type/$selection/82trees_corrected
		echo -en "...with corrected reading frame...and for selection: $selection"
	else
		mafftpath=$type/$selection/60mafft
		alnpath=$type/$selection/70concatenated_exon_alignments
		alnpathselected=$type/$selection/71selected
		treepath=$type/$selection/72trees
		echo -en "...and for selection: $selection"
	fi
fi

#Check necessary file
echo -ne "\nTesting if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
		if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
			echo -e "OK\n"
		else
			echo -e "no alignment files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
			rm -d ../workdir08e 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
else
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
		if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
			echo -e "OK\n"
		else
			echo -e "no alignment files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
			rm -d ../workdir08e 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08g)" ]; then
		echo -e "Directory 'workdir08g' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
fi

#Add necessary scripts and files
#cp $source/readWriteTrees.R .
#Copy list of genes
if [[ $update =~ "yes" ]]; then
	cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
	mv selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
else
	cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
fi

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy/output
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy/input_files
else
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy/output
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy/input_files
fi

#Copy bootrapped gene tree files (if tree=RAxML or tree=FastTree and FastTreeBoot=yes)
if [[ $tree =~ "RAxML" ]]; then
	#Copy RAxML bootstraped trees
	if [[ $update =~ "yes" ]]; then
		cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_${i}_modif${MISSINGPERCENT}.result .
		done
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap* .
	fi
	#Make a list of bootstraped trees
	ls *bootstrap* > bs-files.txt
elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	#Copy FastTree bootstraped trees
	if [[ $update =~ "yes" ]]; then
		cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_${i}_*.boot.fast.trees .
		done
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees .
	fi
	#Make a list of bootstraped trees
	ls *.boot.fast.trees > bs-files.txt
fi

if [[ $cp =~ "yes" ]]; then
	#Removing '_cpDNA' from gene trees in trees.newick
	if [ -d "boot" ]; then
		sed -i.bak 's/_cpDNA//g' *
	fi
fi

# Modify files with bootstrap trees - write as NEXUS and make them compatible with BUCKy
for file in $(cat bs-files.txt); do
	#read newick file and write nexus
	R -e "library(ape); args <- commandArgs(); name <- args[4]; trees<-read.tree(name); write.nexus(trees, file=\"trees.nex\")" $file > /dev/null
	#R --slave -f readWriteTrees.R $file
	#add increasing number after 'TREE'
	perl -pe 's/\bTREE\b/$& . ++$count/ge' trees.nex > trees_modif.nex
	#remove ' * UNTITLED'
	sed -i 's/ \* UNTITLED//' trees_modif.nex
	#remove ' [&U]'
	sed -i 's/ \[\&U\]//' trees_modif.nex
	#change case in 'BEGIN TREES' (to avoid modification in next step)
	sed -i s'/BEGIN TREES/begin trees/' trees_modif.nex
	#change 'TREE' to 'tree rep.'
	sed -i 's/TREE/tree rep\./' trees_modif.nex
	rm trees.nex
	rm $file
	mv trees_modif.nex $file
	#prepare input file for BUCKy
	mbsum -n 0 -o $file.in $file
done

# Run BUCKy
ls *.in > bucky_input.txt
#Select only the last part of the $data (i.e., after the last '/')
datamodif=$(echo $data | awk -F "/" '{ print $NF}')
#a=alpha, n=number of chain steps, k=number of runs, c=number of chains
bucky -i bucky_input.txt -o $datamodif -a $alpha -n $nrbucky -k $nrruns -c $nrchains

#Modify BUCKy output: extract trees in NEXUS format
echo -e "#NEXUS\nbegin trees;" > tree.tre
#Copy 'translate' block everything from 'translate' to ';'
sed -n '/translate/,/;/p' ${datamodif}.concordance >> tree.tre
cat tree.tre > BUCKy_popultree.newick
cat tree.tre > BUCKy_popultreeBL.newick
cat tree.tre > BUCKy_conctree.newick
cat tree.tre > BUCKy_conctreeCF.newick

#Extract 'Population Tree' (the line matching and one following line), change EOL to space and change name to conform NEXUS standards
grep -a1 "^Population Tree:" ${datamodif}.concordance | tr '\n' ' ' | sed 's/ Population Tree:/tree PopulationTree =/' >> BUCKy_popultree.newick
echo -e "\nEND;" >> BUCKy_popultree.newick
#Modify to NEWICK
R -e "library(ape); args <- commandArgs(); name <- args[4]; trees<-read.nexus(name); write.tree(trees, file=\"BUCKy_popultree.tre\")" BUCKy_popultree.newick > /dev/null

#Extract 'Population Tree, With Branch Lengths' (the line matching and one following line), change EOL to space and change name to conform NEXUS standards
grep -a1 "^Population Tree, With Branch Lengths" ${datamodif}.concordance | tr '\n' ' ' | sed 's/ Population Tree, With Branch Lengths In Estimated Coalescent Units:/tree PopulationTreeWithBranchLengthsInEstimatedCoalescentUnits =/' >> BUCKy_popultreeBL.newick
echo -e "\nEND;" >> BUCKy_popultreeBL.newick
#Modify to NEWICK
R -e "library(ape); args <- commandArgs(); name <- args[4]; trees<-read.nexus(name); write.tree(trees, file=\"BUCKy_popultreeBL.tre\")" BUCKy_popultreeBL.newick > /dev/null

#Extract 'Primary Concordance Tree Topology' (the line matching and one following line), change EOL to space and change name to conform NEXUS standards
grep -a1 "^Primary Concordance Tree Topology" ${datamodif}.concordance | tr '\n' ' ' | sed 's/ Primary Concordance Tree Topology:/tree PrimaryConcordanceTreeTopology =/' >> BUCKy_conctree.newick
echo -e "\nEND;" >> BUCKy_conctree.newick
#Modify to NEWICK
R -e "library(ape); args <- commandArgs(); name <- args[4]; trees<-read.nexus(name); write.tree(trees, file=\"BUCKy_conctree.tre\")" BUCKy_conctree.newick > /dev/null

#Extract 'Primary Concordance Tree with' (the line matching and one following line), change EOL to space and change name to conform NEXUS standards
grep -a1 "^Primary Concordance Tree with" ${datamodif}.concordance | tr '\n' ' ' | sed 's/ Primary Concordance Tree with Sample Concordance Factors:/tree PrimaryConcordanceTreewithSampleConcordanceFactors =/' >> BUCKy_conctreeCF.newick
echo -e "\nEND;" >> BUCKy_conctreeCF.newick
#Modify to NEWICK
R -e "library(ape); args <- commandArgs(); name <- args[4]; trees<-read.nexus(name); write.tree(trees, file=\"BUCKy_conctreeCF.tre\")" BUCKy_conctreeCF.newick > /dev/null

rm tree.tre

# Copy results to home
if [[ $update =~ "yes" ]]; then
	cp ${datamodif}* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy/output
	cp BUCKy*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy
	cp *.in $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/BUCKy/input_files
else
	cp ${datamodif}* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy/output
	cp BUCKy*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy
	cp *.in $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/BUCKy/input_files
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08g
fi

echo -e "HybPhyloMaker 8g finished...\n"
