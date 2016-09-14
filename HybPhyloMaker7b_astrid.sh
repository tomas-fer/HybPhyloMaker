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
# *                                   v.1.2.0                                    *
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
	echo -e "\nHybPhyloMaker7b is running on MetaCentrum..."
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
	#module add astrid-1.0
	#module add jdk-7
	module add python-2.7.6-gcc
	module add python27-modules-gcc
	module add python-3.4.1-intel
elif [[ $HOSTNAME == *local* ]]; then
	echo -e "\nHybPhyloMaker7b is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir07b
	cd workdir07b
	#Add necessary modules
	module load bioinformatics/anaconda3/2.3.0
else
	echo -e "\nHybPhyloMaker7b is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir07b
	cd workdir07b
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
	if [ -z "$OUTGROUP" ]; then
		if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick" ]; then
			if [[ $mlbs =~ "yes" ]]; then
				if [[ $tree =~ "RAxML" ]]; then
					if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_*_modif${MISSINGPERCENT}.result 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML'. Exiting..."
							rm -d ../workdir07b/ 2>/dev/null
							exit 3
						fi
					else
						echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
					if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_*_*.boot.fast.trees 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
							rm -d ../workdir07b/ 2>/dev/null
							exit 3
						fi
					else
						echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				fi
			else
				echo -e "OK\n"
			fi
		else
			echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
			rm -d ../workdir07b/ 2>/dev/null
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
							rm -d ../workdir07b/ 2>/dev/null
							exit 3
						fi
					else
						echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
					if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
						if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_*_*.boot.fast.trees 2>/dev/null | wc -w) ]; then
							echo -e "OK\n"
						else
							echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
							rm -d ../workdir07b/ 2>/dev/null
							exit 3
						fi
					else
						echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				fi
			else
				echo -e "OK\n"
			fi
		else
			echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
			rm -d ../workdir07b/ 2>/dev/null
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
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
					if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees 2>/dev/null | wc -w) ]; then
						echo -e "OK\n"
					else
						echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				fi
			else
				echo -e "OK\n"
			fi
		else
			echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
			rm -d ../workdir07b/ 2>/dev/null
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
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
					if [ 0 -lt $(ls $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees 2>/dev/null | wc -w) ]; then
						echo -e "OK\n"
					else
						echo -e "MLBS was requested but no bootstrap replicate trees were found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree'. Exiting..."
						rm -d ../workdir07b/ 2>/dev/null
						exit 3
					fi
				fi
			else
				echo -e "OK\n"
			fi
		else
			echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
			rm -d ../workdir07b/ 2>/dev/null
			exit 3
		fi
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir07b/ 2>/dev/null
		exit 3
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir07b/ 2>/dev/null
		exit 3
	fi
fi
if [ "$(ls -A ../workdir07b)" ]; then
	echo -e "Directory 'workdir07b' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir07b/ 2>/dev/null
	exit 3
fi


#Add necessary programs and files
cp $source/$astridbin .
cp $source/combineboot.py .
cp $source/AMAS.py .

#Copy genetree file
if [[ $update =~ "yes" ]]; then
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick .
		mv trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
	fi
else
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick .
		mv trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
	fi
fi

#Modify labels in gene tree
sed -i.bak 's/XX/-/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
sed -i.bak2 's/YY/_/g' trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick

if [[ $mlbs =~ "yes" ]]; then
	#Copy bootrapped gene tree files (if tree=RAxML or tree=FastTree and FastTreeBoot=yes)
	if [[ $tree =~ "RAxML" ]]; then
		#Make dir for for bootstraped trees
		mkdir boot
		#Copy RAxML bootstraped trees
		if [[ $update =~ "yes" ]]; then
			cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
			for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap*ssembly_${i}_modif${MISSINGPERCENT}.result boot/
			done
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/RAxML_bootstrap* boot/
		fi
		#Make a list of bootstraped trees
		ls boot/*bootstrap* > bs-files
	elif [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
		#Make dir for for bootstraped trees
		mkdir boot
		#Copy FastTree bootstraped trees
		if [[ $update =~ "yes" ]]; then
			cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
			for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt); do
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*ssembly_${i}_*.boot.fast.trees boot/
			done
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree/*.boot.fast.trees boot/
		fi
		#Make a list of bootstraped trees
		ls boot/*.boot.fast.trees > bs-files
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
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid
fi

#Run ASTRID
echo -e "Computing ASTRID multilocus bootstrap..."
#ASTRID -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -m bionj -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
if [[ $tree =~ "RAxML" ]] || [[ $tree =~ "FastTree" && $FastTreeBoot =~ "yes" ]]; then
	if [[ $mlbs =~ "yes" ]]; then
		#Run Astrid bootstrap
		./$astridbin -i trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick -b bs-files -o Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE} > Astrid_boot.log
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
			if [ ! -f $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip ]; then
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
				#Prepare concatenated dataset and transform it to phylip format
				python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip >/dev/null
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
			fi
			if [[ $cp =~ "yes" ]]; then
				#Removing '_cpDNA' from names in alignment
				sed -i.bak 's/_cpDNA//g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
			fi
			mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip concatenated.phylip
			sed -i.bak 's/-/XX/' concatenated.phylip
			sed -i.bak2 's/_/YY/' concatenated.phylip
			#Combine bootstrap with consensus tree
			python ./combineboot.py Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
			mv combinedSupportsTree.tre Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre
			echo
		fi
	else
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
	nw_reroot Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	if [[ $mlbs =~ "yes" ]]; then
		nw_reroot Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
		nw_reroot Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootmajorcons.tre
		if [[ $combine =~ "yes" ]]; then
			nw_reroot Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre $OUTGROUP > tmp && mv tmp Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre
			#Remove possible excessive zeros from combined tree
			sed -i.bak11 's/.0000//g' Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_bootANDcons.tre
		fi
	fi
fi

#Modify labels in Astrid trees
sed -i.bak9 's/-/ /g' Astrid*.tre
sed -i.bak10 's/_/ /g' Astrid*.tre

#Copy results and logs to home
if [[ $update =~ "yes" ]]; then
	cp Astrid*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astrid
else
	cp Astrid*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astrid
fi
echo -e "\nProgress of ASTRID run is written to Astrid.log or Astrid_boot.log (if MLBS was requested)..."

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir07b
fi

echo -e "\nScript HybPhyloMaker7b finished...\n"
