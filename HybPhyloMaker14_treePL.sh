#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker14_treePL
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker14_treePL
#$ -o HybPhyloMaker14_treePL.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                     Script 14 - treePL divergence dating                     *
# *                                   v.1.8.0f                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Requires (in HybSeqSource):
#(1) configuration.txt and
#(2) treepl_wrapper.sh (https://github.com/tongjial/treepl_wrapper)
#Date the species tree (either ExaML, FastTree or Astral4) using penalized likelihood (implemented in treePL, https://github.com/blackrim/treePL)
#collapsed and requisite options not yet implemented!!!

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker14 is running on MetaCentrum..."
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
	module add jdk-1.6.0
	module add treepl-1.0
	#Python (and its modules) are added later in the script
	#module add python-2.7.6-gcc
	#module add python27-modules-gcc
	#module add python-3.4.1-gcc
	module add newick-utils-13042016
	module add r/4.4.0-gcc-10.2.1-oxdi5pz
	#module add debian11/compat #necessary for R-3.4.3
	#module add R-3.4.3-gcc
	#module add debian8-compat
	#module add p4 #do not load before running 'python3 AMAS.py'
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker14 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir14
	cd workdir14
	#Add necessary modules
	module load java/1.7
	module load bioinformatics/anaconda3/5.1 #python3 and NewickUtilities
	module load tools/R/3.4.1
	#module load bioinformatics/newickutilities/0.0
	#module load bioinformatics/p4/ #???
else
	echo -e "\nHybPhyloMaker14 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	export LD_LIBRARY_PATH=/usr/local/lib64
	#Make and enter work directory
	mkdir -p workdir14
	cd workdir14
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

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir14)" ]; then
		echo -e "Directory 'workdir14' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir14 2>/dev/null
		exit 3
	fi
fi

#Add necessary programs and files
if [ -f "$source/configuration.txt" ]; then
	cp $source/configuration.txt .
else
	echo -e "The file 'configuration.txt' is missing in HybSeqSource. Exiting...\n"
	exit 3
fi
if [ -f "$source/treepl_wrapper.sh" ]; then
	cp $source/treepl_wrapper.sh .
else
	echo -e "The script 'treepl_wrapper.sh' is missing in HybSeqSource. Exiting...\n"
	exit 3
fi

#Write log
logname=HPM14
echo -e "HybPhyloMaker14: treePL divergence dating" > ${logname}.log
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

#Set folder
if [[ $tpltree =~ "ExaML" ]]; then
	tpltf=concatenatedExaML
elif [[ $tpltree =~ "Astral4" ]]; then
	tpltf=Astral4
elif [[ $tpltree =~ "FastTree" ]]; then
	tpltf=concatenated
fi

#Check if the results exist with other species trees
if [[ $update =~ "yes" ]]; then
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/treePL/bstrees" ]; then
		tpltb=concatenatedExaML
		echo -e "...treePL results for bootstrap trees probably found associated with '$tpltb' tree in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/treePL/bstrees' folder."
		echo -e "Check it and proceed with running the script 14b (available treePL results will be used to annotate dated $tpltree species tree)."
		echo -e "Exiting...\n"
		exit 3
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral4/treePL/bstrees" ]; then
		tpltb=Astral4
		echo -e "...treePL results for bootstrap trees probably found associated with '$tpltb' tree in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral4/treePL/bstrees' folder."
		echo -e "Check it and proceed with running the script 14b (available treePL results will be used to annotate dated $tpltree species tree)."
		echo -e "Exiting...\n"
		exit 3
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/treePL/bstrees" ]; then
		tpltb=concatenated
		echo -e "...treePL results for bootstrap trees probably found associated with '$tpltb' tree in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/treePL/bstrees' folder."
		echo -e "Check it and proceed with running the script 14b (available treePL results will be used to annotate dated $tpltree species tree)."
		echo -e "Exiting...\n"
		exit 3
	else
		echo -e "No treePL results for bootstrap trees found. Continue running the script...\n"
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/treePL/bstrees" ]; then
		tpltb=concatenatedExaML
		echo -e "...treePL results for bootstrap trees probably found associated with '$tpltb' tree in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/treePL/bstrees' folder."
		echo -e "Check it and proceed with running the script 14b (available treePL results will be used to annotate dated $tpltree species tree)."
		echo -e "Exiting...\n"
		exit 3
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral4/treePL/bstrees" ]; then
		tpltb=Astral4
		echo -e "...treePL results for bootstrap trees probably found associated with '$tpltb' tree in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral4/treePL/bstrees' folder."
		echo -e "Check it and proceed with running the script 14b (available treePL results will be used to annotate dated $tpltree species tree)."
		echo -e "Exiting...\n"
		exit 3
	elif [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/treePL/bstrees" ]; then
		tpltb=concatenated
		echo -e "...treePL results for bootstrap trees probably found associated with '$tpltb' tree in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/treePL/bstrees' folder."
		echo -e "Check it and proceed with running the script 14b (available treePL results will be used to annotate dated $tpltree species tree)."
		echo -e "Exiting...\n"
		exit 3
	else
		echo -e "No treePL results for bootstrap trees found. Continue running the script...\n"
	fi
fi

#Create folder for results & Copy species tree
echo -e "Creating folder for results...\n"
if [[ $tpltree =~ "ExaML" ]]; then
	tplt=ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	if [[ $update =~ "yes" ]]; then
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/treePL
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	else
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/treePL
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	fi
elif [[ $tpltree =~ "Astral4" ]]; then
	tplt=Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	if [[ $update =~ "yes" ]]; then
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral4/treePL
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral4/Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	else
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral4/treePL
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral4/Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
	fi
elif [[ $tpltree =~ "FastTree" ]]; then
	tplt=concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
	if [[ $update =~ "yes" ]]; then
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/treePL
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre .
	else
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/treePL
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre .
	fi
else
	echo -e "'tpltree=' has to be one of ExaML, Astral4 or FastTree. Exiting...\n"
	exit 3
fi

if [[ $tpltree =~ "Astral4" ]]; then
	#Modify Astral4 tree (replace ' ' back to '-' and '_')
	sed -i 's/ \([^ ]*\) / \1_/g' Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre #replace every second occurrence of ' ' by '_'
	sed -i 's/ /-/g' Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre #replace all spaces by '-'
elif [[ $tpltree =~ "FastTree" ]]; then
	#Modify concatenated FastTree tree (replace ' ' back to '-' and '_')
	sed -i 's/ \([^ ]*\) / \1_/g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre #replace every second occurrence of ' ' by '_'
	sed -i 's/ /-/g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre #replace all spaces by '-'
fi

#Root tree
nw_reroot -s ${tplt} $OUTGROUP > tmp && mv tmp ${tplt}

#Preapare for running treePL on bootstrap trees
if [[ $tplbs =~ "yes" ]]; then
	#Check if treePL bootstrap analysis already exists
	if [[ $update =~ "yes" ]]; then
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees" ]; then
			echo -e "The folder '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees' already exists. Delete it or rename before running this script again. Exiting..."
			exit 3
		fi
	else
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees" ]; then
			echo -e "The folder '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees' already exists. Delete it or rename before running this script again. Exiting..."
			exit 3
		fi
	fi
	#Copy gene tree file & root trees with OUTGROUP
	echo -e "Copying gene trees & separating & creating jobs...\n"
	if [[ $update =~ "yes" ]]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	fi
	#Only gene trees with OUTGROUP
	grep ${OUTGROUP} ${treefile} > tmp && mv tmp ${treefile}
	nw_reroot -s ${treefile} $OUTGROUP > bstrees
	rm ${treefile}
	#Split trees to single files (each line to separate file)
	awk '{filename = sprintf("tree%d.tre", NR); print >filename; close(filename)}' bstrees
	#Create individual folders, copy input files there
	if [[ $update =~ "yes" ]]; then
		mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees
	else
		mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees
	fi
	for i in $(ls tree*.tre | cut -d'.' -f1); do
		if [[ $update =~ "yes" ]]; then
			mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees/${i}
			cp ${i}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees/${i}
		else
			mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees/${i}
			cp ${i}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees/${i}
		fi
	done
	#Prepare job files
	touch $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees/submitTreePLjobs.sh
	for group in $(ls tree*.tre | cut -d'.' -f1); do
		echo '#!/bin/bash' >> ${group}.sh
		echo '#----------------MetaCentrum----------------' >> ${group}.sh
		echo '#PBS -l walltime=24:00:00' >> ${group}.sh
		echo '#PBS -l select=1:ncpus=4:mem=1gb:scratch_local=1gb' >> ${group}.sh
		echo '#PBS -j oe' >> ${group}.sh
		echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> ${group}.sh
		echo '#PBS -N treePL_for_'"${group}" >> ${group}.sh
		echo '#Add necessary modules' >> ${group}.sh
		echo 'module add treepl-1.0' >> ${group}.sh
		echo 'cd $SCRATCHDIR' >> ${group}.sh
		echo 'cp '"$source"'/treepl_wrapper.sh .' >> ${group}.sh
		echo 'cp '"$source"'/configuration.txt .' >> ${group}.sh
		echo 'chmod +x treepl_wrapper.sh' >> ${group}.sh
		echo 'cp '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}"'/'"${group}"'.tre .' >> ${group}.sh
		echo '#Run PLtree within treepl_wrapper' >> ${group}.sh
		echo './treepl_wrapper.sh configuration.txt '"${group}"'.tre '"${group}"' > treePL_'"${group}"'.log' >> ${group}.sh
		echo '#Copy results to home' >> ${group}.sh
		echo 'cp configure* '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}" >> ${group}.sh
		echo 'cp cv* '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}" >> ${group}.sh
		echo 'cp out_dates* '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}" >> ${group}.sh
		echo 'cp prime* '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}" >> ${group}.sh
		echo 'cp treepl*.tre* '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}" >> ${group}.sh
		echo 'cp *.log '"$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}"'/update/species_trees/'"${modif}${tpltf}"'/treePL/bstrees/'"${group}" >> ${group}.sh
		chmod +x ${group}.sh
		if [[ $update =~ "yes" ]]; then
			cp ${group}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees/
			echo 'qsub '"${group}"'.sh' >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees/submitTreePLjobs.sh
		else
			cp ${group}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees/
			echo 'qsub '"${group}"'.sh' >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL/bstrees/submitTreePLjobs.sh
		fi
	done
	#Make the submitter executable
	chmod +x $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees/submitTreePLjobs.sh
fi

#Run main treePL analysis
if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/${tpltree}_treePLresult.tre" ]; then
	echo -e "The file '${tpltree}_treePLresult.tre' already exist in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/'. Delete it or rename before running this script again. No treePL analysis now...\n"
else
	echo -e "Running treePL for the main species tree (${tpltree})...\n"
	chmod 755 treepl_wrapper.sh
	./treepl_wrapper.sh configuration.txt ${tplt} treePLresult > treePL_${tpltree}.log
	#Modify names
	mv treepl_treePLresult.tre ${tpltree}_treePLresult.tre
	mv treepl_treePLresult.tre.r8s ${tpltree}_treePLresult.tre.r8s
	mv out_dates.tre.r8s ${tpltree}_dates.tre.r8s
	mv out_dates.tre ${tpltree}_dates.tre
	#Delete some files
	rm ${tplt} treepl_wrapper.sh
	# Copy results to home
	if [[ $update =~ "yes" ]]; then
		cp * $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL
		cp treePL_${tpltree}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL
	else
		cp * $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL
		cp treePL_${tpltree}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpltf}/treePL
	fi
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
	rm -r workdir14
fi

if [[ $tplbs =~ "yes" ]]; then
	echo -e "Go to '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpltf}/treePL/bstrees' and run 'submitTreePLjobs.sh'\n"
fi
echo -e "\nScript HybPhyloMaker14 finished...\n"
