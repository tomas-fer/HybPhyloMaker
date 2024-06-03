#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=24gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8j_superQ_network
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8j_superQ_network
#$ -o HybPhyloMaker8j_superQ_network.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                    Script 08j - SuperQ network in Spectre                    *
# *                                   v.1.8.0d                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2024 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute SuperQ network from selected gene trees using Spectre
#Take gene trees specified in trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick in 'species_trees' folder

#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6b_FastTree_for_selected.sh or HybPhyloMaker6a_RAxML_for_selected.sh to create gene trees
#(3) HybPhyloMaker7_roottrees.sh to create the file containing all gene trees with BS values removed

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8j is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add spectre-1.1.15
	module add debian9-compat
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8j is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08j
	cd workdir08j
	#Add necessary modules
else
	echo -e "\nHybPhyloMaker8j is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08j
	cd workdir08j
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
		rm -d ../workdir08j 2>/dev/null
		exit 3
	fi
else
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
		echo -e "OK\n"
	else
		echo -e "no gene trees file found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/'. Run first HybPhyloMaker7_roottrees.sh. Exiting..."
		rm -d ../workdir08j 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08j 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08j 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08j)" ]; then
		echo -e "Directory 'workdir08j' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08j 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8j
echo -e "HybPhyloMaker8j: SuperQ network in Spectre" > ${logname}.log
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
for set in data selection cp corrected update tree; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
fi


# Copy gene trees
echo -e "Copying gene trees...\n"
if [[ $update =~ "yes" ]]; then
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
else
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
fi

#Compute SuperQ network
echo -e "Computing SuperQ network...\n"
if [[ ! $location == "1" ]]; then
	#copy SuperQ binary dir
	cp -r $source/spectre-1.1.5 .
	#export JAVA_OPTS='-Xmx80g'
	spectre-1.1.5/bin/superq -o output.net -y JOptimizer -b balanced trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick > spectre.log
else
	export JAVA_OPTS='-Xmx80g'
	superq -o output.net -y JOptimizer -b balanced trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick > spectre.log
fi
mv output.net SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex

#Removing '_cpDNA' from names in network
sed -i.bak 's/_cpDNA//g' SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex

#Renaming sample names
sed -i.bak 's/XX/-/g' SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex
sed -i.bak 's/YY/_/g' SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex

#Plot network using phangorn in R
module add r/4.4.0-gcc-10.2.1-ssuwpvb
export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"
cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex net.nex
R -q -e "library(phangorn);nnet<-read.nexus.splits('net.nex');a=8;pdf('net.pdf');par(mar=c(a,a,a,a),xpd=T);plot(as.networx(nnet),cex=0.3,tip.color='firebrick4',edge.width=0.3,edge.lty=1,direction='axial');dev.off();write.nexus.networx(as.networx(nnet),file='n.nex')" >> R.log 2>&1
mv n.nex SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network.nex

#Crop white margins in PDF using GhostScript
module add ghostscript
gs -o null -sDEVICE=bbox net.pdf 2>out #reports crop box coordinates
#take the four numbers and add/subtract a value ('add')
add=10
crop1=$(grep HiRes out | cut -d' ' -f2) #left margin
crop1x=$(echo "$crop1 - $add" | bc)
crop2=$(grep HiRes out | cut -d' ' -f3) #bottom margin
crop2x=$(echo "$crop2 - $add" | bc)
crop3=$(grep HiRes out | cut -d' ' -f4) #right margin
crop3x=$(echo "$crop3 + $add" | bc)
crop4=$(grep HiRes out | cut -d' ' -f5) #upper margin
crop4x=$(echo "$crop4 + $add" | bc)
crop=$(echo $crop1x $crop2x $crop3x $crop4x)
gs -o netcrop.pdf -sDEVICE=pdfwrite -dAutoRotatePages=/None -dUseCropBox=true -c "[/CropBox [$crop] /PAGES pdfmark" -f net.pdf
mv net.pdf SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network.pdf
mv netcrop.pdf SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network_cropped.pdf

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network.nex $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network_cropped.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
	cp spectre.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
	cp R.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
else
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_splits.nex $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network.nex $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
	cp SuperQNetwork_${MISSINGPERCENT}_${SPECIESPRESENCE}_network_cropped.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
	cp spectre.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
	cp R.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SuperQ
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SuperQ
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08j
fi

echo -e "\nHybPhyloMaker8j finished...\n"
