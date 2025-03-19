#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker14b_treePL_BSsummary
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker14b_treePL_BSsummary
#$ -o HybPhyloMaker14_treePL_BSsummary.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *               Script 14b - summary of treePL on bootstrap trees              *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker14b is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add jdk/8
	module add beast2
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker14b is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir14b
	cd workdir14b
	#Add necessary modules
	module load java/1.8
	module load beast2?
else
	echo -e "\nHybPhyloMaker14b is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir14b
	cd workdir14b
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

#Set folder
if [[ $tpltree =~ "ExaML" ]]; then
	tpltf=concatenatedExaML
elif [[ $tpltree =~ "Astral4" ]]; then
	tpltf=Astral4
elif [[ $tpltree =~ "FastTree" ]]; then
	tpltf=concatenated
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir14b)" ]; then
		echo -e "Directory 'workdir14b' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir14b 2>/dev/null
		exit 3
	fi
fi

#Add necessary programs and files
if [ -f "$source/configuration.txt" ]; then
	cp $source/configuration.txt .
else
	echo -e "The file 'configuration.txt' is missing in HybSeqSource. Exiting...\n"
	rm -d ../workdir14b 2>/dev/null
	exit 3
fi

#Write log
logname=HPM14b
echo -e "HybPhyloMaker14b: summary of treePL on bootstrap trees" > ${logname}.log
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
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree FastTreeBoot OUTGROUP collapse requisite tpltree; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

echo -e "\nSettings for divergence dating:" >> ${logname}.log
cat configuration.txt >> ${logname}.log

if [[ $requisite =~ "yes" ]]; then
	echo -e "\nList of requisite samples" >> ${logname}.log
	echo $requisitetaxa | tr '|' '\n' >> ${logname}.log
fi
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Check if the results exist with other species trees
if [[ $update =~ "yes" ]]; then
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/treePL/bstrees" ]; then
		tpltb=concatenatedExaML
		echo -e "...treePL results for bootstrap trees found in '$tpltb' folder...\n"
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral4/treePL/bstrees" ]; then
		tpltb=Astral4
		echo -e "...treePL results for bootstrap trees found in '$tpltb' folder...\n"
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/treePL/bstrees" ]; then
		tpltb=concatenated
		echo -e "...treePL results for bootstrap trees found in '$tpltb' folder...\n"
	else
		echo -e "No treePL results for bootstrap trees found. Run first script 14. Exiting...\n"
		rm -d ../workdir14b 2>/dev/null
		exit 3
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/treePL/bstrees" ]; then
		tpltb=concatenatedExaML
		echo -e "...treePL results for bootstrap trees found in '$tpltb' folder...\n"
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral4/treePL/bstrees" ]; then
		tpltb=Astral4
		echo -e "...treePL results for bootstrap trees found in '$tpltb' folder...\n"
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/treePL/bstrees" ]; then
		tpltb=concatenated
		echo -e "...treePL results for bootstrap trees found in '$tpltb' folder...\n"
	else
		echo -e "No treePL results for bootstrap trees found. Run first script 14. Exiting...\n"
		rm -d ../workdir14b 2>/dev/null
		exit 3
	fi
fi

#Copy results of treePL on bootstrap trees
echo -e "Copying treePL results on bootstrap trees...\n"
if [[ $update =~ "yes" ]]; then
	cp -r $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltb}/treePL/bstrees/ .
else
	cp -r $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltb}/treePL/bstrees/ .
fi

#Copy main treePL tree
echo -e "Copying treePL result on $tpltree species tree...\n"
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/${tpltree}_treePLresult.tre .
else
	cp * $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/${tpltree}_treePLresult.tre .
fi

#Create a file with results of all bootstrap trees
echo -e "Creating a file with results of all bootstrap trees...\n"
for i in $(ls -d bstrees/tree* | grep -v sh | cut -d'/' -f2); do
	cat bstrees/${i}/treepl_${i}.tre >> ${tpltree}_bstreepl.tre
done

#Summarize dated trees and create a tree with intervals
echo -e "Treeannotator: summarizing dated bootstrap trees and creating a tree with intervals...\n"
treeannotator -heights median -target ${tpltree}_treePLresult.tre ${tpltree}_bstreepl.tre ${tpltree}_treePL_BSintervalsMedian.tre 2> treeanotator.log

#Copy results back to home
if [[ $update =~ "yes" ]]; then
	cp ${tpltree}_bstreepl.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL
	cp ${tpltree}_treePL_BSintervalsMedian.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL
	cp treeanotator.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL
else
	cp ${tpltree}_bstreepl.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL
	cp ${tpltree}_treePL_BSintervalsMedian.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL
	cp treeanotator.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL
fi

#Finish log & copy home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir14b
fi

echo -e "\nScript HybPhyloMaker14b finished...\n"

