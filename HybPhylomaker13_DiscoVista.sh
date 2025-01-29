#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=4gb
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
# *                                   v.1.8.0c                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
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
	#. /packages/run/modules-2.0/init/bash
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

#Write log
logname=HPM13
echo -e "HybPhyloMaker13: DiscoVista" > ${logname}.log
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
for set in data selection cp corrected OUTGROUP; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi


#clone DiscoVista repository
git clone https://github.com/esayyari/DiscoVista
export WS_HOME=`pwd`

#copy list of folders to compare
cp $source/comparelist.txt .

#method 1 (species tree)
mkdir species parameters
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
	dirn=${str}.${ty}.${mis}.${pres}.${corr}.${sl}.${tr}.${up}.${col}-FNA
	mkdir species/$dirn
	#check for Astral
	if [ ! -z $(grep Astral <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/Astral_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/Astral_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for Astrid
	if [ ! -z $(grep Astrid <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/Astrid_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/Astrid_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for MRL
	if [ ! -z $(grep MRL <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/MRL_${mis}_${pres}.tre species/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/MRL_${mis}_${pres}_${col}.tre species/$dirn/estimated_species_tree.tree
		fi
		sed -i 's/ /_/g' species/$dirn/estimated_species_tree.tree #replace all ' ' occurrences by '-'
	fi
	#check for concatenated
	if [ ! -z $(grep "concatenated$" <<< $i) ]; then
		if [[ $col =~ "nocol" ]]; then
			cp ${path}${i}/concatenated${mis}_${pres}.fast.tre species/$dirn/estimated_species_tree.tree
		else
			cp ${path}${i}/concatenated${mis}_${pres}_${col}.fast.tre species/$dirn/estimated_species_tree.tree
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
done

#Run DiscoVista
$WS_HOME/DiscoVista/src/utils/discoVista.py -c parameters/clade-defs.txt -p species/ -t 0.95 -m 0 -o results

#Make dir for results
mkdir -p $path/DiscoVista/speciestree
#Copy results back to home
cp -r species/ $path/DiscoVista/speciestree
cp -r parameters/ $path/DiscoVista/speciestree
cp -r results/ $path/DiscoVista/speciestree
cp comparelist.txt $path/DiscoVista/speciestree

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/DiscoVista

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
