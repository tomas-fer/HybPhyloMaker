#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus=10:mem=16gb:scratch_local=2gb
#PBS -j oe
#PBS -N HybPhyloMaker8l_SNaQ
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 10
#$ -q mThC.q
#$ -l mres=2G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8l_SNaQ
#$ -o HybPhyloMaker8l_SNaQ.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                  Script 08l - SNaQ network in PhyloNetworks                  *
# *                                   v.1.8.0a                                   *
# *                           Roman Ufimov & Tomas Fer                           *
# * Dept. of Botany, Charles University, Prague, Czech Republic, 2023            *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute SNaQ network from selected gene trees using PhyloNetworks
#Take gene trees specified in trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick in 'species_trees' folder
#Take Astral species tree

#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6b_FastTree_for_selected.sh or HybPhyloMaker6a_RAxML_for_selected.sh to create gene trees
#(3) HybPhyloMaker7_roottrees.sh to create the file containing all gene trees with BS values removed
#(4) HybPhyloMaker8a_astral.sh to create Astral species tree

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8l is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add r/4.0.0-gcc
	module add julia
	module add newick-utils-13042016
	module add python36-modules-gcc #to add cairosvg
	#module add debian9-compat
	#set nr processors
	cpu=$TORQUE_RESC_PROC
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8l is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08j
	cd workdir08j
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #adds NewickUtilities
	module load tools/R/3.4.1
	#module load bioinformatics/newickutilities/0.0
	#julia
	#cairosvg
else
	echo -e "\nHybPhyloMaker8l is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08l
	cd workdir08l
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
echo -ne "\n\nTesting if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
		echo -e "OK\n"
	else
		echo -e "no gene trees file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/'. Run first HybPhyloMaker7_roottrees.sh. Exiting..."
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
else
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
		echo -e "OK\n"
	else
		echo -e "no gene trees file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/'. Run first HybPhyloMaker7_roottrees.sh. Exiting..."
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNaQ" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNaQ' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNaQ" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNaQ' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08l)" ]; then
		echo -e "Directory 'workdir08l' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8l
echo -e "HybPhyloMaker8l: SNaQ network in PhyloNetworks" > ${logname}.log
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "run on MetaCentrum: $PBS_O_HOST" >> ${logname}.log
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
for set in data selection cp corrected update tree hstart hmax; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
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

#Define Astral tree name
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

#Remove '.tre' from species tree name
sptree=$(cut -d'.' -f1 <<< $astraltree)

#Copy species tree
echo -e "\nCopying trees..."
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/${astraltree} .
else
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/${astraltree} .
fi

#Copy gene trees (file with all gene trees)
if [[ $update =~ "yes" ]]; then
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	fi
else
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	fi
fi

#Modify Astral tree (replace ' ' back to '-' and '_')
sed -i 's/ \([^ ]*\) / \1_/g' ${sptree}.tre #replace every second occurrence of ' ' by '_'
sed -i 's/ /-/g' ${sptree}.tre #replace all spaces by '-'
#Add a space at the end of the line
sed -i 's/$/ /' ${sptree}.tre
#Removing '_cpDNA' from names
sed -i.bak 's/_cpDNA//g' ${sptree}.tre

#Reroot Astral tree with $OUTGROUP
nw_reroot -s ${sptree}.tre $OUTGROUP > tmp && mv tmp ${sptree}.tre

#Modify labels in gene tree
sed -i.bak 's/XX/-/g' $treefile
sed -i.bak2 's/YY/_/g' $treefile
#Removing '_cpDNA' from names
sed -i.bak 's/_cpDNA//g' $treefile

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNaQ
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNaQ
fi

#Rename input files
mv $treefile trees.newick
mv ${sptree}.tre sptree.tre

#Compute SNaQ network
echo -e "Computing SNaQ network...\n"
cp $source/runSNaQ.jl .
julia runSNaQ.jl "$OUTGROUP" $hstart $hmax $cpu > runSNaQ.log

#Rename 
mv Rplots.pdf NetworkScoresPlot.pdf
#Create PDF from all SVG
echo -e "Transforming SVG to PDF...\n"
for i in $(ls *.svg | cut -d'.' -f1); do
	cairosvg ${i}.svg -o ${i}.pdf
done

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp *.{log,out,networks,csv,pdf,svg,txt} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNaQ
else
	cp *.{log,out,networks,csv,pdf,svg,txt} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNaQ
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNaQ
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNaQ
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08l
fi

echo -e "\nHybPhyloMaker8l finished...\n"
