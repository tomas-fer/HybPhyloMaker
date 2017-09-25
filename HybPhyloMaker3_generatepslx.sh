#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=6:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker3_generate_pslx
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker3_generate_pslx
#$ -o HybPhyloMaker3_generate_pslx.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *         Script 03 - Process consensus after mapping, make pslx files         *
# *                                   v.1.5.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2017 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

# Input: Consensus sequences from Geneious: must be named consensus.fasta or consensus_cpDNA.fasta 
# it is multiple fasta file with names,
# e.g., Camptandra-latifolia_S4-all-no-dups_assembled_to_Curcuma_exons_reference_400Ns__consensus_sequence)
# prepared in /storage/$server/home/$LOGNAME/data/30consensus/

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker3 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add blat-suite-34
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker3 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir -p workdir03
	cd workdir03
	#Add necessary modules
	module load bioinformatics/blat/36x1
else
	echo -e "\nHybPhyloMaker3 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir -p workdir03
	cd workdir03
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
else
	echo -e "Working with exons\n"
	type="exons"
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [[ $cp =~ "yes" ]]; then
	if [ -f "$path/$type/30consensus/consensus_cpDNA.fasta" ]; then
		if [ -f "$source/$cpDNACDS" ]; then
			echo -e "OK\n"
		else
			echo -e "'$cpDNACDS' is missing in 'HybSeqSource'. Exiting...\n"
			rm -d ../workdir03/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/$type/30consensus/consensus_cpDNA.fasta' is missing. Exiting...\n"
		rm -d ../workdir03/ 2>/dev/null
		exit 3
	fi
else
	if [ -f "$path/$type/30consensus/consensus.fasta" ]; then
		if [ -f "$source/$probes" ]; then
			echo -e "OK\n"
		else
			echo -e "'$probes' is missing in 'HybSeqSource'. Exiting...\n"
			rm -d ../workdir03/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/$type/30consensus/consensus.fasta' is missing. Exiting...\n"
		rm -d ../workdir03/ 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [ -d "$path/$type/40contigs" ]; then
	echo -e "Directory '$path/$type/40contigs' already exists. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir03/ 2>/dev/null
	exit 3
else
	if [ -d "$path/$type/50pslx" ]; then
		echo -e "Directory '$path/$type/50pslx' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir03/ 2>/dev/null
		exit 3
	else
		if [[ ! $location == "1" ]]; then
			if [ "$(ls -A ../workdir03)" ]; then
				echo -e "Directory 'workdir03' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir03/ 2>/dev/null
				exit 3
			fi
		fi
	fi
fi

#Copy fasta from home folder to scratch/workdir
cp -r $path/$type/30consensus/* .
#Make a new folder for results
mkdir -p $path/$type
mkdir $path/$type/40contigs

#-----------------------GENEIOUS CONSENSUS SEQUENCE MODIFICATION-----------------------
echo -en "Parsing consensus sequence..."
#Modify Windows EOLs to Unix EOLs (i.e., LF only)
if [[ $cp =~ "yes" ]]; then
	sed -i.bak 's/\x0D$//' consensus_cpDNA.fasta
else
	sed -i.bak 's/\x0D$//' consensus.fasta
fi
#Remove trailing '?'s, remove unwanted part of the file name (everything after the second '-'), add '_consensus_sequence' to all headers
#and split multiple fasta from Geneious into individual fasta sequences
if [[ $cp =~ "yes" ]]; then
	cat consensus_cpDNA.fasta | sed 's/^?*//' | cut -d"-" -f1,2 | sed '/>/s/.*/&_consensus_cpDNAsequence/' | awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'
	rm consensus_cpDNA.fasta
	rm -f consensus.fasta
else
	cat consensus.fasta | sed 's/^?*//' | cut -d"-" -f1,2 | sed '/>/s/.*/&_consensus_sequence/' | awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'
	rm consensus.fasta
	rm -f consensus_cpDNA.fasta
fi
#Former version
#cat consensus.fasta | sed 's/^?*//' | sed "s/-all-no-dups_assembled_to_${pseudoref}_//" | awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'

#Make a file with names of a all fasta files in a folder. Name is modified to consist only of first two parts of the name separated by '_' 
ls *.fasta | cut -d'_' -f 1,2 > listOfConsensusFiles.txt
#A loop for preparing assemblies (consensus sequences) of individual exons from Geneious consensus
if [[ $cp =~ "yes" ]]; then
	for file in $(cat listOfConsensusFiles.txt)
	do
		#Modify Windows EOLs to Unix EOLs (i.e., LF only)
		#sed -i 's/.$//' $file\_consensus_cpDNAsequence.fasta
		
		#Delete first line in a fasta file, i.e. header
		sed -i.bak 1d $file\_consensus_cpDNAsequence.fasta
		
		#Replace sequence of '?' by new line (\n), put bash variable ($file) to 'val' which is available to print with awk,
		#print '>Contig'+number(NR;increased by one each step)+species name (val), then to next line print sequence
		cat $file\_consensus_cpDNAsequence.fasta | tr -s '?' '\n' | awk -v val=$file '{ print ">Contig" NR "_" val "\n" $1 m}' > $file\_contigs_cpDNA.fas
		#Copy data from scratch to home folder
		cp $file\_contigs_cpDNA.fas $path/$type/40contigs
	done
else
	for file in $(cat listOfConsensusFiles.txt)
	do
		#Modify Windows EOLs to Unix EOLs (i.e., LF only)
		#sed -i 's/.$//' $file\_consensus_sequence.fasta
		
		#Delete first line in a fasta file, i.e. header
		sed -i.bak 1d $file\_consensus_sequence.fasta
		
		#Replace sequence of '?' by new line (\n), put bash variable ($file) to 'val' which is available to print with awk,
		#print '>Contig'+number(NR;increased by one each step)+species name (val), then to next line print sequence
		cat $file\_consensus_sequence.fasta | tr -s '?' '\n' | awk -v val=$file '{ print ">Contig" NR "_" val "\n" $1 m}' > $file\_contigs.fas
		#Copy data from scratch to home folder
		cp $file\_contigs.fas $path/$type/40contigs
	done
fi
echo -e "finished\n"

#-----------------------BLAT ASSEMBLIES TO REFERENCE-----------------------
echo -e "Generating pslx files using BLAT...\n"
#Copy other transcriptome/genome data from home to scratch/workdir (must be named with suffix *.fas)
if [ "$othersource" != "" ] && [ "$othersource" != "NO" ]; then 
	cp -r $othersourcepath/* .
fi

if [[ $cp =~ "yes" ]]; then
	#Copy cpDNA reference
	cp -r $source/$cpDNACDS .
else
	#Copy reference
	cp -r $source/$probes .
fi
#Make a new folder for results
mkdir $path/$type/50pslx

#Make a list of all files with contigs
ls *.fas > contig_names.txt

#A loop to process all contig files specified in contig_names.txt
for contigfile in $(cat contig_names.txt)
do
	echo -e "\nProcessing $contigfile..."
	if [[ $cp =~ "yes" ]]; then
		blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $cpDNACDS $contigfile ${contigfile}.pslx
	else
		blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $probes $contigfile ${contigfile}.pslx
	fi
	cp $contigfile.pslx $path/$type/50pslx
done

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir03
fi

echo -e "\nScript HybPhyloMaker3 finished...\n"
