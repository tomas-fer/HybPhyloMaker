#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=8gb
#PBS -N HybPhyloMaker7b_Astrid
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=2G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker7b_Astrid
#$ -o HybPhyloMaker7b_Astrid.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                       Script 07b - Astrid species tree                       *
# *                                   v.1.1.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using Astrid methods from trees saved in single gene tree file (with *.newick suffix)
#Take trees from 72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Run first
#(1) HybPhyloMaker4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker5a_RAxML_for_selected.sh or HybPhyloMaker5b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPhyloMaker6_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values


#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
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
elif [[ $HOSTNAME == *local* ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07b
	cd workdir07b
	#Add necessary modules
	
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07b
	cd workdir07b
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

echo -e "\nScript HybPhyloMaker7b is running..."

#Copy genetree file
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
else
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
fi

if [[ $mlbs =~ "yes" ]]; then
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
fi
#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid
else
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid
fi

#Run ASTRID
echo -e "\nComputing Astrid tree..."
$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
#ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -m bionj -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	if [[ $mlbs =~ "yes" ]]; then
		#Run Astrid bootstrap
		$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre
		#Extract last tree (main species tree based on non-bootstrapped dataset + bootstrap values mapped on it)
		cat Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre | tail -n1 > Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
	fi
fi

#Modify labels in Astrid trees
sed -i 's/XX/-/g' Astrid*.tre
sed -i 's/YY/_/g' Astrid*.tre

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp Astrid*.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid
else
	cp Astrid*.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid
fi
#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07b
fi

echo -e "\nScript HybPhyloMaker7b finished..."
