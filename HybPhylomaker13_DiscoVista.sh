#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker13_DiscoVista
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker13_DiscoVista
#$ -o HybPhyloMaker13_DiscoVista.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                            Script 13 - DiscoVista                            *
# *                                   v.1.8.1                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker13 is running on MetaCentrum..."
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
	module add R-3.4.3-gcc
	module add python27-modules-gcc #adds also biopython and DendroPy
	#module add debian8-compat
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker13 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir13
	cd workdir13
	#Add necessary modules
	module load tools/R/3.4.1
else
	echo -e "\nHybPhyloMaker13 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir13
	cd workdir13
fi

#clone DiscoVista repository
git clone https://github.com/esayyari/DiscoVista
export WS_HOME=`pwd`

# R
# install.packages(c("reshape","reshape2","ggplot2","plyr","scales","ape","optparse"))
# pip install dendropy #Python2

#copy list of folders to compare
cp $source/comparelist.txt .

#method 0, 1, 5  (species tree)
mkdir species parameters genetrees rel_frequency
#copy annotations
cp $source/annotation.txt parameters/
#create clade-defs.txt file
$WS_HOME/DiscoVista/src/utils/generate_clade-defs.py parameters/annotation.txt parameters/clade-defs.txt

#copy species trees defined in 'comparelist.txt' to respective folders
for i in $(cat comparelist.txt); do
	if [ ! -z $(grep exon <<< $i) ]; then ty=exon; fi
	if [ ! -z $(grep cp <<< $i) ]; then ty=cp; fi
	if [ ! -z $(grep corrected <<< $i) ]; then corr=corr; else corr=orig; fi
	if [ ! -z $(grep RAxML <<< $i) ]; then tr=RAxML; fi
	if [ ! -z $(grep FastTree <<< $i) ]; then tr=FastTree; fi
	if [ ! -z $(grep update <<< $i) ]; then up=update; else up=all; fi
	if [ ! -z $(grep collapsed <<< $i) ]; then col=$(tr '/' '\n' <<< $i | grep collapsed); else col=nocol; fi
	if [ $(tr '/' '\n' <<< $i | sed "s/_corrected/_corrected\n/" | sed 's/species_trees//' | tr '_' '\n' | sed '/^$/d' | grep trees -n | cut -d ':' -f 1) -eq 2 ]; then sl=nosel; fi
	if [ $(tr '/' '\n' <<< $i | sed "s/_corrected/_corrected\n/" | sed 's/species_trees//' | tr '_' '\n' | sed '/^$/d' | grep trees -n | cut -d ':' -f 1) -eq 3 ]; then sl=$(tr '/' '\n' <<< $i | sed '/^$/d' | sed -n 2p); fi
	mis=$(tr '/' '\n' <<< $i | grep [0-9]tree | sed 's/82trees_corrected//' | sed 's/72trees//' | cut -d '_' -f 1)
	pres=$(tr '/' '\n' <<< $i | grep [0-9]tree | sed 's/82trees_corrected//' | sed 's/72trees//' | cut -d '_' -f 2)
	if [ ! -z $(grep Astral <<< $i) ]; then str=Astral; fi
	if [ ! -z $(grep Astrid <<< $i) ]; then str=Astrid; fi
	if [ ! -z $(grep MRL <<< $i) ]; then str=MRL; fi
	if [ ! -z $(grep "concatenated$" <<< $i) ]; then str=concatFT; fi
	if [ ! -z $(grep concatenatedExaML <<< $i) ]; then str=concatExaML; fi
	dirn=${str}.${ty}.${mis}.${pres}.${corr}.${sl}.${tr}.${up}.${col}
	mkdir species/$dirn
	mkdir genetrees/$dirn
	mkdir rel_frequency/$dirn

	# fix i to not contain SPECIENSTREE FOLDER
	speciestree_folder=$(sed 's![^/]*$!!' <<< ${i})

	cp ${path}${speciestree_folder}/trees.newick genetrees/$dirn/genetrees.tree
	sed -i 's/-/_/g' genetrees/$dirn/genetrees.tree #replace all '-' occurrences by '-'

	cp ${path}${speciestree_folder}/trees.newick rel_frequency/$dirn/estimated_gene_trees.tree
	sed -i 's/-/_/g' rel_frequency/$dirn/estimated_gene_trees.tree #replace all '-' occurrences by '-'

	#check for Astral
	if [ ! -z $(grep Astral <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/Astral_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/Astral_${mis}_${pres}.tre rel_frequency/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/Astral_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/Astral_${mis}_${pres}_${col}.tre rel_frequency/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for Astrid
	if [ ! -z $(grep Astrid <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/Astrid_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/Astrid_${mis}_${pres}.tre rel_frequency/$dirn/estimated_species_tree.tree		
		else
			cp ${path}${i}/Astrid_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/Astrid_${mis}_${pres}_${col}.tre rel_frequency/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for MRL
	if [ ! -z $(grep MRL <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/MRL_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/MRL_${mis}_${pres}.tre rel_frequency/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/MRL_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/MRL_${mis}_${pres}_${col}.tre rel_frequency/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for concatenated
	if [ ! -z $(grep "concatenated$" <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/concatenated${mis}_${pres}.fast.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/concatenated${mis}_${pres}.fast.tre rel_frequency/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/concatenated${mis}_${pres}_${col}.fast.tre species/$dirn/estimated_species_tree.tree
			cp ${path}${i}/concatenated${mis}_${pres}_${col}.fast.tre rel_frequency/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for concatenatedExaML
	if [ ! -z $(grep concatenatedExaML <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/ExaML_bootstrap_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/ExaML_bootstrap_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
		fi
	fi
	sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	sed -i 's/ /_/g' rel_frequency/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'

	
	# prepare data for m1
	# write each gene tree to single file
	mkdir -p genetrees/input
	num=1
	while read line
	do
		mkdir -p genetrees/input/tree_${num}
		mkdir genetrees/input/tree_${num}/${dirn}
	    echo "$line" > genetrees/input/tree_${num}/${dirn}/estimated_gene_trees.tree
	    let num=num+1
	done < genetrees/$dirn/genetrees.tree

	# run m5
	mkdir rel_frequency/$dirn/results

	#Run DiscoVista M5
	$WS_HOME/DiscoVista/src/utils/discoVista.py -c parameters/clade-defs.txt -p rel_frequency/$dirn/ -t 90 -m 5 -o rel_frequency/$dirn/results -a parameters/annotation.txt -g Outgroup
done

#Run DiscoVista M0
$WS_HOME/DiscoVista/src/utils/discoVista.py -c parameters/clade-defs.txt -p species/ -t 0.90 -m 0 -o results
#Run DiscoVista M1
$WS_HOME/DiscoVista/src/utils/discoVista.py -c parameters/clade-defs.txt -p genetrees/input/ -t 90 -m 1 -o genetrees/results -a parameters/annotation.txt -k 20

#Make dir for results
mkdir -p $path/DiscoVista/speciestree
#Copy results back to home
cp -r species/ $path/DiscoVista/speciestree
cp -r parameters/ $path/DiscoVista/speciestree
cp -r results/ $path/DiscoVista/speciestree
cp comparelist.txt $path/DiscoVista/speciestree

p -r genetrees/ $path/DiscoVista/
cp comparelist.txt $path/DiscoVista/genetrees
cp -r parameters/ $path/DiscoVista/genetrees

cp -r rel_frequency/ $path/DiscoVista/
cp -r parameters/ $path/DiscoVista/rel_frequency
cp comparelist.txt $path/DiscoVista/rel_frequency

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir13
fi

echo -e "\nScript HybPhyloMaker13 finished...\n"
