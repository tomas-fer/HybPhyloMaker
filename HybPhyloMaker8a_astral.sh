#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker8a_Astral
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8a_Astral
#$ -o HybPhyloMaker8a_Astral.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                       Script 08a - Astral species tree                       *
# *                                   v.1.6.6                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using ASTRAL methods from trees saved in a single gene tree file (with *.newick suffix)
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6a_RAxML_for_selected.sh or HybPhyloMaker6b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPhyloMaker7_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#Works also for trees after update, requisite taxa selection and collapsing (see HybPhyloMaker9_update_trees.sh and HybPhyloMaker10_requisite_collapse.sh)
#Calculation of multilocus bootstrap does not work for trees after collapsing (but works also for trees with requisite samples only)

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8a is running on MetaCentrum..."
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
	#Python (and its modules) are added later in the script
	#module add python-2.7.6-gcc
	#module add python27-modules-gcc
	module add python-3.4.1-gcc
	module add newick-utils-13042016
	module add debian8-compat
	#module add p4 #do not load before running 'python3 AMAS.py'
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8a is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08a
	cd workdir08a
	#Add necessary modules
	module load java/1.7
	module load bioinformatics/anaconda3/5.1 #python3 and NewickUtilities
	#module load bioinformatics/newickutilities/0.0
	#module load bioinformatics/p4/ #???
else
	echo -e "\nHybPhyloMaker8a is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08a
	cd workdir08a
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
	echo -en "...with corrected reading frame"
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
fi

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection"
else
	echo -e ""
fi

if [[ ! $collapse -eq "0" ]]; then
	echo -e "...and with trees with branches below ${collapse} BS collapsed"
else
	if [[ $requisite =~ "no" ]]; then
		echo -e "\n"
	fi
fi

if [[ $requisite =~ "yes" ]]; then
	echo -e "...and only with trees with requisite taxa present\n"
else
	echo -e "\n"
fi

#Settings for collapsed and requisite selection
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		modif=with_requisite/collapsed${collapse}/
		treefile=trees_with_requisite_collapsed${collapse}.newick
	else
		modif=with_requisite/
		if [ -z "$OUTGROUP" ]; then
			treefile=trees_with_requisite.newick
		else
			treefile=trees_rooted_with_requisite.newick
		fi
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		modif=collapsed${collapse}/
		treefile=trees_collapsed${collapse}.newick
	else
		modif=""
		if [ -z "$OUTGROUP" ]; then
			treefile=trees.newick
		else
			treefile=trees_rooted.newick
		fi
	fi
fi

#Check necessary file
if [[ $requisite =~ "no" ]] && [[ $collapse -eq "0" ]];then
	echo -ne "Testing if input data are available..."
	if [[ $update =~ "yes" ]]; then
		if [ -z "$OUTGROUP" ]; then
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick" ]; then
				if [[ $mlbs =~ "yes" ]]; then
					if [[ $tree =~ "RAxML" ]]; then
						if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
							if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_*_modif${MISSINGPERCENT}.result 2>/dev/null | wc -w) ]; then
								echo -e "OK\n"
							else
								echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML'. Exiting..."
								rm -d ../workdir08a 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
							if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_*_*.boot.fast.trees 2>/dev/null | wc -w) ]; then
								echo -e "OK\n"
							else
								echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
								rm -d ../workdir08a 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08a 2>/dev/null
				exit 3
			fi
		else
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
				if [[ $mlbs =~ "yes" ]]; then
					if [[ $tree =~ "RAxML" ]]; then
						if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
							if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_*_modif${MISSINGPERCENT}.result 2>/dev/null | wc -w) ]; then
								echo -e "OK\n"
							else
								echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML'. Exiting..."
								rm -d ../workdir08a 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
							if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_*_*.boot.fast.trees 2>/dev/null | wc -w) ]; then
								echo -e "OK\n"
							else
								echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
								rm -d ../workdir08a 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08a 2>/dev/null
				exit 3
			fi
		fi
	else
		if [ -z "$OUTGROUP" ]; then
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick" ]; then
				if [[ $mlbs =~ "yes" ]]; then
					if [[ $tree =~ "RAxML" ]]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap* 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML'. Exiting..."
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08a 2>/dev/null
				exit 3
			fi
		else
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
				if [[ $mlbs =~ "yes" ]]; then
					if [[ $tree =~ "RAxML" ]]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap* 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML'. Exiting..."
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
							rm -d ../workdir08a 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08a 2>/dev/null
				exit 3
			fi
		fi
	fi
	
	#Test if folder for results exits
	if [[ $update =~ "yes" ]]; then
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral" ]; then
			echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir08a 2>/dev/null
			exit 3
		fi
	else
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral" ]; then
			echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir08a 2>/dev/null
			exit 3
		fi
	fi
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08a)" ]; then
		echo -e "Directory 'workdir08a' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08a 2>/dev/null
		exit 3
	fi
fi

#Add necessary programs and files
cp $source/$astraljar .
cp -r $source/lib .
cp $source/combineboot.py .
cp $source/AMAS.py .

#Copy genetree file
if [[ $update =~ "yes" ]]; then
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
		#mv trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	fi
else
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
		#mv trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	fi
fi

#Modify labels in gene tree
sed -i.bak 's/XX/-/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
sed -i.bak2 's/YY/_/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick

if [[ $collapse -eq "0" ]];then
	if [[ $mlbs =~ "yes" ]]; then
		#Copy bootrapped gene tree files (if tree=RAxML or tree=FastTree and FastTreeBoot=yes)
		if [[ $tree =~ "RAxML" ]]; then
			#Make dir for for bootstraped trees
			mkdir boot
			#Copy RAxML bootstraped trees
			if [[ $update =~ "yes" ]]; then
				if [[ $requisite =~ "yes" ]]; then
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/selected_genes_with_requisite.txt .
					mv selected_genes_with_requisite.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt
				else
					cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
				fi
				for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_${i}_modif${MISSINGPERCENT}.result boot/
				done
			else
				if [[ $requisite =~ "yes" ]]; then
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/selected_genes_with_requisite.txt .
					for i in $(cat selected_genes_with_requisite.txt); do
						cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_${i}_modif${MISSINGPERCENT}.result boot/
					done
				else
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap* boot/
				fi
			fi
			#Make a list of bootstrapped trees
			ls boot/*bootstrap* > bs-files
		elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
			#Make dir for for bootstrapped trees
			mkdir boot
			#Copy FastTree bootstrapped trees
			if [[ $update =~ "yes" ]]; then
				if [[ $requisite =~ "yes" ]]; then
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/selected_genes_with_requisite.txt .
					mv selected_genes_with_requisite.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt
				else
					cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
				fi
				for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_${i}_*.boot.fast.trees boot/
				done
			else
				if [[ $requisite =~ "yes" ]]; then
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/selected_genes_with_requisite.txt .
					for i in $(cat selected_genes_with_requisite.txt); do
						cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_${i}_*.boot.fast.trees boot/
					done
				else
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees boot/
				fi
			fi
			#Make a list of bootstrapped trees
			ls boot/*.boot.fast.trees > bs-files
		fi
	fi
fi

if [[ $cp =~ "yes" ]]; then
	#Removing '_cpDNA' from gene trees in trees.newick
	if [ -d "boot" ]; then
		sed -i.bak 's/_cpDNA//g' boot/*
	fi
fi

#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral
fi

#Run ASTRAL
echo -e "Computing ASTRAL tree..."
if [[ $location == "1" ]]; then
	java -jar $astraljar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre 2> Astral.log
elif [[ $location == "2" ]]; then
	java -d64 -server -XX:MaxHeapSize=4g -jar $astraljar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre 2> Astral.log
else
	java -jar $astraljar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre 2> Astral.log
fi

#(Re)root a final Astral species tree with $OUTGROUP
if [ -n "$OUTGROUP" ]; then
	nw_reroot -s Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $OUTGROUP > tmp && mv tmp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	echo
fi
#Make a copy of the main Astral tree for future combination(s)
cp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}main.tre
sed -i.bak 's/-/XX/g' Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}main.tre
sed -i.bak2 's/_/YY/g' Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}main.tre

#Modify labels in Astral trees
sed -i.bak 's/-/ /g' Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
sed -i.bak2 's/_/ /g' Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre

#Rename Astral trees
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite_collapsed${collapse}.tre
	else
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.tre
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_collapsed${collapse}.tre
	else
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	fi
fi

#Copy species tree and log to home
if [[ $update =~ "yes" ]]; then
	cp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/${astraltree}
	cp Astral.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral
else
	cp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/${astraltree}
	cp Astral.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral
fi

if [[ $collapse -eq "0" ]];then
	if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
		if [[ $mlbs =~ "yes" ]]; then
			#Run Astral bootstrap
			echo -e "\nComputing ASTRAL multilocus bootrap..."
			if [[ $location == "1" ]]; then
				java -jar $astraljar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre 2> Astral_boot.log
			elif [[ $location == "2" ]]; then
				java -d64 -server -XX:MaxHeapSize=4g -jar $astraljar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre 2> Astral_boot.log
			else
				java -jar $astraljar -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre 2> Astral_boot.log
			fi
			#Extract last tree (main species tree based on non-bootstrapped dataset + bootstrap values mapped on it)
			cat Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre | tail -n1 > Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
			#Extract last but one tree (majority rule consensus tree of bootstrap replicate trees)
			cat Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre | tail -n2 | head -n1 > Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
			#Leave only bootstrap trees in *_allbootstraptrees.tree (delete last two trees extracted above)
			tac Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre | sed '1,2d' | tac > tmp && mv tmp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre
			#(Re)root a final Astral species tree with $OUTGROUP
			if [ -n "$OUTGROUP" ]; then
				nw_reroot -s Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre $OUTGROUP > tmp && mv tmp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
				nw_reroot -s Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre $OUTGROUP > tmp && mv tmp Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
				echo
			fi
			sed -i.bak3 's/-/XX/g' Astral*.tre
			sed -i.bak4 's/_/YY/g' Astral*.tre
			#Make combined trees
			if [[ $combine =~ "yes" ]]; then
				echo -e "\nCombining support values from main, bootstrap and bootstrap consensus trees to one tree..."
				#Copy alignment (for sample names) and rename it
				if [ ! -f $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip ]; then
					#Copy list of genes
					if [[ $update =~ "yes" ]]; then
						cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
						mv selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
					else
						cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
					fi
					#Copy and modify selected FASTA files
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
					#Add python3 module if on MetaCentrum
					# if [[ $PBS_O_HOST == *".cz" ]]; then
						# module add python-3.4.1-intel
					# fi
					#Prepare concatenated dataset and transform it to phylip format
					python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip >/dev/null
				else
					cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
				fi
				if [[ $cp =~ "yes" ]]; then
					#Removing '_cpDNA' from names in alignment
					sed -i.bak 's/_cpDNA//g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
				fi
				mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip concatenated.phylip
				sed -i.bak 's/-/XX/' concatenated.phylip
				sed -i.bak2 's/_/YY/' concatenated.phylip
				#Remove python3 module and load p4 module if on MetaCentrum
				if [[ $PBS_O_HOST == *".cz" ]]; then
					module rm python-3.4.1-gcc
					#module add python-2.7.6-gcc
					#module add python27-modules-gcc
					module add p4
				elif [[ $HOSTNAME == compute-*-*.local ]]; then
					module unload bioinformatics/anaconda3
					#module load bioinformatics/p4?
				fi
				#Combine basic Astral tree with bootstrap tree
				if [[ $PBS_O_HOST == *".cz" ]]; then
					#on Metacentrum, p4 within 'combineboot.py' requires python 2
					python ./combineboot.py Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}main.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
				else
					#locally, p4 within 'combineboot.py' requires python 3
					python3 ./combineboot.py Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}main.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
				fi
				mv combinedSupportsTree.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_mainANDboot.tre
				#Combine with bootstrap consensus tree
				if [[ $PBS_O_HOST == *".cz" ]]; then
					#on Metacentrum, p4 within 'combineboot.py' requires python 2
					python ./combineboot.py Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_mainANDboot.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
				else
					#locally, p4 within 'combineboot.py' requires python 3
					python3 ./combineboot.py Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_mainANDboot.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
				fi
				mv combinedSupportsTree.tre Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_mainANDbootANDcons.tre
				echo
			fi
			#Modify labels in Astral trees
			sed -i.bak5 's/XX/ /g' Astral*.tre
			sed -i.bak6 's/YY/ /g' Astral*.tre
			#Copy results and log to home
			if [[ $update =~ "yes" ]]; then
				cp Astral*boot*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral
				cp Astral_boot.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral
			else
				cp Astral*boot*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral
				cp Astral_boot.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral
			fi
		fi
	fi
fi
echo -e "\nProgress of ASTRAL run is written to Astral.log and Astral_boot.log (if MLBS was requested)..."

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08a
fi

echo -e "\nScript HybPhyloMaker8a finished...\n"
