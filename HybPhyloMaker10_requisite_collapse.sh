#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker10_requisite_collapse
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker10_requisite_collapse
#$ -o HybPhyloMaker10_requisite_collapse.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *  Script 10 - Select trees with requisite taxa, collapse unsupported branches *
# *                                   v.1.8.0b                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Select only gene trees and selected loci containing 'requisite' samples
#If requested, also collapse unsupported branches in gene trees (according to 'collapse' in settings.cfg) using TreeCollapseCL4.jar
#Produces
# - a NEWICK file only with trees containing requisite samples (trees_with_requisite.newick)
# - a text file with alignment file names containing requisite samples (selected_genes_with_requisite.txt)
# - a NEWICK file with trees with unsupported branches collapsed (trees_with_requisite_collapsed${collapse}.newick)

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker10_requisite_collapse is running on MetaCentrum..."
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
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker10_requisite_collapse is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir10
	cd workdir10
	#Add necessary modules
	module load java/1.7
else
	echo -e "\nHybPhyloMaker10_requisite_collapse is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir10
	cd workdir10
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
		echo -e "...with corrected reading frame"
	else
		mafftpath=$type/60mafft
		alnpath=$type/70concatenated_exon_alignments
		alnpathselected=$type/71selected
		treepath=$type/72trees
		echo -e ""
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/$selection/61mafft_corrected
		alnpath=$type/$selection/80concatenated_exon_alignments_corrected
		alnpathselected=$type/$selection/81selected_corrected
		treepath=$type/$selection/82trees_corrected
		echo -e "...with corrected reading frame...and for selection: $selection"
	else
		mafftpath=$type/$selection/60mafft
		alnpath=$type/$selection/70concatenated_exon_alignments
		alnpathselected=$type/$selection/71selected
		treepath=$type/$selection/72trees
		echo -e "...and for selection: $selection"
	fi
fi

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection"
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

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir10)" ]; then
		echo -e "Directory 'workdir10' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir10 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM10
echo -e "HybPhyloMaker10: select trees with requisite taxa, collapse unsupported branches" > ${logname}.log
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
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree OUTGROUP collapse requisite requisitetaxa; do
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
cp $source/TreeCollapseCL4.jar .

if [[ $requisite =~ "no" ]] && [[ $collapse -eq "0" ]]; then
	echo "Nothing to do..."
fi

#Copy genetree file
if [[ $requisite =~ "yes" ]]; then
	if [[ $update =~ "yes" ]]; then
		if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/trees_rooted_with_requisite.newick" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/trees_rooted_with_requisite.newick .
		else
			if [ -z "$OUTGROUP" ]; then
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
				mv trees.newick trees_rooted.newick
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees_rooted.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
			fi
		fi
	else
		if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/trees_rooted_with_requisite.newick" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/trees_rooted_with_requisite.newick .
		else
			if [ -z "$OUTGROUP" ]; then
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
				mv trees.newick trees_rooted.newick
			else
				cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees_rooted.newick .
				mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
			fi
		fi
	fi
else
	if [[ $update =~ "yes" ]]; then
		if [ -z "$OUTGROUP" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees.newick .
			mv trees.newick trees_rooted.newick
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees_rooted.newick .
		fi
	else
		if [ -z "$OUTGROUP" ]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees.newick .
			mv trees.newick trees_rooted.newick
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees_rooted.newick .
		fi
	fi
fi

#Select only trees with requisite taxa and copy the file to home
if [[ $requisite =~ "yes" ]]; then
	if [ ! -f "trees_rooted_with_requisite.newick" ]; then
		echo -e "Selecting trees with requisite taxa..."
		grep -E "$requisitetaxa" trees_rooted.newick > trees_rooted_with_requisite.newick
		#Copy selected trees to home
		if [[ $update =~ "yes" ]]; then
			if [ -z "$OUTGROUP" ]; then
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/trees_with_requisite.newick
			else
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
			fi
		else
			if [ -z "$OUTGROUP" ]; then
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/trees_with_requisite.newick
			else
				cp trees_rooted_with_requisite.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
			fi
		fi
	fi
fi

#Make a list of genes containing requisite taxa in alignment
if [[ $requisite =~ "yes" ]]; then
	echo -e "Preparing list of genes containing requisite taxa..."
	#grep only filenames of alignments (option '-l') containing requisite taxa (sed command extract only filenames without extension from the full path)
	grep -El "$requisitetaxa" $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*modif* | sed -r "s/.+\/(.+)\..+/\1/" > genes_with_requisite_ALL.txt
	if [[ $update =~ "yes" ]]; then
		#Copy list of selected genes
		cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		#Add '_' before and after each numbers
		awk '{ print "_" $0 "_" }' selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt > tmp && mv tmp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt
		#Take only selected genes with requisite taxa
		grep -F -f selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt < genes_with_requisite_ALL.txt > selected_genes_with_requisite.txt
		cut -d'_' -f2 selected_genes_with_requisite.txt > tmp && mv tmp selected_genes_with_requisite.txt
		cp selected_genes_with_requisite.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
	else
		#Copy list of selected genes
		cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
		#Take only numbers of selected genes
		cut -d'_' -f2 selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt > tmp && mv tmp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
		#Add '_' before and after each numbers
		awk '{ print "_" $0 "_" }' selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt > tmp && mv tmp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
		#Take only selected genes with requisite taxa
		grep -F -f selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt < genes_with_requisite_ALL.txt > selected_genes_with_requisite.txt
		cut -d'_' -f2 selected_genes_with_requisite.txt > tmp && mv tmp selected_genes_with_requisite.txt
		cp selected_genes_with_requisite.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
	fi
fi

#Collapse trees using 'TreeCollapseCL4'
if [[ ! $collapse -eq "0" ]]; then
	echo -e "Collapsing tree branches with BS below ${collapse}..."
	unset i
	mkdir trees
	#Put all trees one-by-one to folder 'tree'
	if [[ $requisite =~ "yes" ]]; then
		for x in $(<trees_rooted_with_requisite.newick); do
			echo "$x" > trees/tree$((++i)).tre
		done
		if [[ $location == "2" ]]; then
			java -d64 -server -XX:MaxHeapSize=1g -jar TreeCollapseCL4.jar -b ${collapse} -d trees/ 1>/dev/null
		else
			java -jar TreeCollapseCL4.jar -b ${collapse} -d trees/ 1>/dev/null
		fi
		cat trees/*coll* > trees_with_requisite_collapsed${collapse}.newick
		if [[ $update =~ "yes" ]]; then
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/collapsed${collapse}
			cp trees_with_requisite_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/collapsed${collapse}
		else
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/collapsed${collapse}
			cp trees_with_requisite_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/collapsed${collapse}
		fi
	else
		for x in $(<trees_rooted.newick); do
			echo "$x" > trees/tree$((++i)).tre
		done
		if [[ $location == "2" ]]; then
			java -d64 -server -XX:MaxHeapSize=1g -jar TreeCollapseCL4.jar -b ${collapse} -d trees/ 1>/dev/null
		else
			java -jar TreeCollapseCL4.jar -b ${collapse} -d trees/ 1>/dev/null
		fi
		cat trees/*coll* > trees_collapsed${collapse}.newick
		if [[ $update =~ "yes" ]]; then
			mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/collapsed${collapse}
			cp trees_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/collapsed${collapse}
		else
			mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/collapsed${collapse}
			cp trees_collapsed${collapse}.newick $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/collapsed${collapse}
		fi
	fi
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ ! $collapse -eq "0" ]]; then
	if [[ $requisite =~ "yes" ]]; then
		if [[ $update =~ "yes" ]]; then
			cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite/collapsed${collapse}
		else
			cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite/collapsed${collapse}
		fi
	else
		if [[ $update =~ "yes" ]]; then
			cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/collapsed${collapse}
		else
			cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/collapsed${collapse}
		fi
	fi
else
	if [[ $requisite =~ "yes" ]]; then
		if [[ $update =~ "yes" ]]; then
			cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/with_requisite
		else
			cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/with_requisite
		fi
	fi
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir10
fi

echo -e "\nScript HybPhyloMaker10_requisite_collapse finished...\n"
