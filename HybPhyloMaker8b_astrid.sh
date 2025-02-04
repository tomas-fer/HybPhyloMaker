#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker8b_Astrid
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=2G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8b_Astrid
#$ -o HybPhyloMaker8b_Astrid.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                       Script 08b - Astrid species tree                       *
# *                                   v.1.8.0g                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using Astrid methods from gene trees trees saved in a single gene tree file (with *.newick suffix)
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6a_RAxML_for_selected.sh or HybPhyloMaker6b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPhyloMaker7_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#Works also for trees after update, requisite taxa selection and collapsing (see HybPhyloMaker9_update_trees.sh and HybPhyloMaker10_requisite_collapse.sh)
#Calculation of multilocus bootstrap does not work for trees after collapsing (but works also for trees with requisite samples only)

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8b is running on MetaCentrum..."
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	#module add astrid-1.0
	#module add jdk-7
	#module add python-2.7.6-gcc
	#module add python27-modules-gcc
	#module add python-3.4.1-gcc
	module add newick-utils-13042016
	#module add debian8-compat
	#module add p4 #do not load before running 'python3 AMAS.py'
elif [[ $HOSTNAME == *local* ]]; then
	echo -e "\nHybPhyloMaker8b is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08b
	cd workdir08b
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #python3 and NewickUtilities
	module load bioinformatics/astrid/1.4
	#module load bioinformatics/p4/ #???
else
	echo -e "\nHybPhyloMaker8b is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08b
	cd workdir08b
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

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection"
else
	echo -e ""
fi

if [[ ! $collapse -eq "0" ]]; then
	echo -e "...and with trees with branches below ${collapse} BS collapsed"
else
	if [[ $requisite =~ "no" ]]; then
		echo -e ""
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
								rm -d ../workdir08b/ 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
							if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_*_*.boot.fast.trees 2>/dev/null | wc -w) ]; then
								echo -e "OK\n"
							else
								echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
								rm -d ../workdir08b/ 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08b/ 2>/dev/null
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
								rm -d ../workdir08b/ 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
							if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_*_*.boot.fast.trees 2>/dev/null | wc -w) ]; then
								echo -e "OK\n"
							else
								echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
								rm -d ../workdir08b/ 2>/dev/null
								exit 3
							fi
						else
							echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08b/ 2>/dev/null
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
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08b/ 2>/dev/null
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
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
							rm -d ../workdir08b/ 2>/dev/null
							exit 3
						fi
					fi
				else
					echo -e "OK\n"
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08b/ 2>/dev/null
				exit 3
			fi
		fi
	fi
	#Test if folder for results exits
	if [[ $update =~ "yes" ]]; then
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid" ]; then
			echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir08b/ 2>/dev/null
			exit 3
		fi
	else
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid" ]; then
			echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir08b/ 2>/dev/null
			exit 3
		fi
	fi
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08b)" ]; then
		echo -e "Directory 'workdir08b' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08b/ 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8b
echo -e "HybPhyloMaker8b: ASTRID species tree" > ${logname}.log
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "Job run on MetaCentrum: $PBS_JOBID" >> ${logname}.log
	echo -e "From: $PBS_O_HOST" >> ${logname}.log
	echo -e "Host: $HOSTNAME" >> ${logname}.log
	echo -e "$PBS_NUM_NODES node(s) with $PBS_NCPUS core(s)" >> ${logname}.log
	memM=$(bc <<< "scale=2; $(echo $PBS_RESC_MEM) / 1024 / 1024 ")
	memG=$(bc <<< "scale=2; $(echo $PBS_RESC_MEM) / 1024 / 1024 / 1024 ")
	if (( $(echo $memG 1 | awk '{if ($1 < $2) print 1;}') )); then
		echo -e "Memory: $memM Mb" >> ${logname}.log
	else
		echo -e "Memory: $memG Gb" >> ${logname}.log
	fi
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "run on Hydra: $HOSTNAME" >> ${logname}.log
else
	echo -e "local run: "`hostname`"/"`whoami` >> ${logname}.log
fi
echo -e "\nBegin:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
echo -e "\nSettings" >> ${logname}.log
if [[ $PBS_O_HOST == *".cz" ]]; then
	printf "%-25s %s\n" `echo -e "\nServer:\t$server"` >> ${logname}.log
fi
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree FastTreeBoot OUTGROUP collapse requisite mlbs combine; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [[ $requisite =~ "yes" ]]; then
	echo -e "\nList of requisite samples" >> ${logname}.log
	echo $requisitetaxa | tr '|' '\n' >> ${logname}.log
fi
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Add necessary programs and files
cp $source/$astridbin .
chmod +x $astridbin
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

#Check if OUTGROUP does exist in gene tree files
if [ ! -z "$OUTGROUP" ]; then #only if outgroup is set
	if [ -z $(grep $OUTGROUP trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick | head -n1) ]; then
		echo -e "$OUTGROUP was not found in gene tree file. Exiting..."
		exit 3
	fi
fi

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
		sed -i.bak3 's/_cpDNA//g' boot/*
	fi
fi

#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astrid
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astrid
fi

#Run ASTRID
#ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -m bionj -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
if [[ $collapse -eq "0" ]];then
	if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
		if [[ $mlbs =~ "yes" ]]; then
			#Run Astrid bootstrap
			echo -e "Computing ASTRID multilocus bootstrap..."
			if [[ $location == "2" ]]; then
				ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE} > Astrid_boot.log
			else
				./$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE} > Astrid_boot.log
			fi
			#Remove "'" from resulting trees
			sed -i.bak4 "s/'//g" Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}*
			#Rename resulting trees
			mv Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE} Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
			mv Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.bs_tree Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
			mv Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.bs_consensus Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
			mv Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.bs Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre
			#Make combined trees
			if [[ $combine =~ "yes" ]]; then
				echo -e "\nCombining support values from main, bootstrap and bootstrap consensus trees to one tree..."
				sed -i.bak5 's/-/XX/g' Astrid*.tre
				sed -i.bak6 's/_/YY/g' Astrid*.tre
				
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
				#Remove python3 module and load p4 module if on MetaCentrum or Hydra
				if [[ $PBS_O_HOST == *".cz" ]]; then
					#module rm python-3.4.1-gcc
					#module add python-2.7.6-gcc
					#module add python27-modules-gcc
					#module add debian8-compat
					module add py-p4phylogenetics/20240606
				elif [[ $HOSTNAME == compute-*-*.local ]]; then
					module unload bioinformatics/anaconda3
					#module load bioinformatics/p4?
				fi
				#Combine bootstrap with consensus tree
				if [[ $PBS_O_HOST == *".cz" ]]; then
					#on Metacentrum, p4 within 'combineboot.py' requires python 2
					python ./combineboot.py Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
					module rm debian8-compat
				else
					#locally, p4 within 'combineboot.py' requires python 3
					python3 ./combineboot.py Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
				fi
				mv combinedSupportsTree.tre Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre
				echo
			fi
		else
			if [[ $location == "2" ]]; then
				ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE} > Astrid_boot.log
			else
				./$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre > Astrid.log
			fi
			#Remove "'" from resulting trees
			sed -i.bak "s/'//g" Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
		fi
	else
		echo -e "Computing ASTRID tree..."
		./$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre > Astrid.log
		#Remove "'" from resulting trees
		sed -i.bak "s/'//g" Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	fi
else
	echo -e "Computing ASTRID tree..."
	./$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre > Astrid.log
	#Remove "'" from resulting trees
	sed -i.bak "s/'//g" Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
fi

#Modify labels in Astrid trees
sed -i.bak7 's/XX/-/g' Astrid*.tre
sed -i.bak8 's/YY/_/g' Astrid*.tre

#(Re)root a final Astrid species tree with $OUTGROUP
if [ -n "$OUTGROUP" ]; then
	nw_reroot -s Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	if [[ $collapse -eq "0" ]];then
		if [[ $mlbs =~ "yes" ]]; then
			nw_reroot -s Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
			nw_reroot -s Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
			if [[ $combine =~ "yes" ]]; then
				nw_reroot -s Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre
				#Remove possible excessive zeros from combined tree
				sed -i.bak11 's/.0000//g' Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre
			fi
		fi
	fi
fi

#Rename Astrid trees
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		astridtree=Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite_collapsed${collapse}.tre
	else
		astridtree=Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.tre
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		astridtree=Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_collapsed${collapse}.tre
	else
		astridtree=Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	fi
fi

#Modify labels in Astrid trees
sed -i.bak9 's/-/ /g' Astrid*.tre
sed -i.bak10 's/_/ /g' Astrid*.tre

#Copy results and logs to home
if [[ $collapse -eq "0" ]];then
	if [[ $update =~ "yes" ]]; then
		cp Astrid*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astrid
		cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astrid
	else
		cp Astrid*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astrid
		cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astrid
	fi
else
	if [[ $update =~ "yes" ]]; then
		cp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astrid/${astridtree}
		cp Astrid.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astrid
	else
		cp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astrid/${astridtree}
		cp Astrid.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astrid
	fi
fi
echo -e "Progress of ASTRID run is written to Astrid.log or Astrid_boot.log (if MLBS was requested)..."

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astrid
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astrid
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08b
fi

echo -e "\nScript HybPhyloMaker8b finished...\n"
