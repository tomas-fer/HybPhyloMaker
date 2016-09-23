#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=4:ppn=6
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=8gb
#PBS -N HybPhyloMaker0a_datadownloadprepare
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 24
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker0a_datadownloadprepare
#$ -o HybPhyloMaker0a_datadownloadprepare.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                     Script 00a - Download & prepare data                     *
# *                                   v.1.2.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


# Download fastq files from BaseSpace (optional), rename them and move to sample folders
# i.e., prepare raw reads for running HybPhyloMaker
# Works for PE reads (i.e., 2 fastq files per sample) in the following format:
# samplename_sampleID_L001_R1_001.fastq.gz
# samplename_sampleID_L001_R2_001.fastq.gz
#
# Needs two files (token_header.txt, renamelist.txt) to be in the home folder (at desired data server if running on cluster)
# Download from BaseSpace is parallelized

#####################################################################################################################
# Usage:                                                                                                            #
# 1. You must obtain 'token' from Illumina BaseSpace (how to do this see steps 1-5 at                               #
# https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader                        #
# 2. Save this token to text file (token_header.txt) with one line text:                                            #
# header = "x-access-token: <your-token-here>                                                                       #
# 3. Login to BaseSpace via web browser and get IDs for                                                             #
# - forward read (R1) of the first sample in a run                                                                  #
# - reverse read (R2) of the last sample in a run                                                                   #
# How to do this: go to (via clicking) Projects -> <project-name> -> Samples -> <sample-name> -> <file>.fastq.gz    #
# Look at the address which should looks like                                                                       #
# https://basespace.illumina.com/sample/28555179/files/tree/Z001_S1_L001_R1_001.fastq.gz?id=2016978377              #
# desired ID is the last number                                                                                     #
# 4. prepare file renamelist.txt with two columns (desired sample name and 'Sample_Name' from BaseSpace)            #
#   (Sample_Name is first part of file name at BaseSpace), e.g.                                                     #
# genus1-species1_S001	Z001                                                                                        #
# genus1-species2_S002	Z002                                                                                        #
# etc.                                                                                                              #
# To ensure smooth processing with subsequent scripts in HybPhyloMaker the file names must follow the convention    #
# genus-species_code                                                                                                #
#                                                                                                                   #
#####################################################################################################################

if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker0a is running on MetaCentrum...\n"
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
	module add parallel
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker0a is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir00_dataprep
	cd workdir00_dataprep
	#Add necessary modules
	module load tools/gnuparallel/20160422
else
	echo -e "\nHybPhyloMaker0a is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir00_dataprep
	cd workdir00_dataprep
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -d "$path/10rawreads" ]; then
	echo -e "Directory '$path/10rawreads' already exists. Exiting...\n"
	rm -d ../workdir00_dataprep/ 2>/dev/null
	exit 3
else
	if [ -f "../renamelist.txt" ]; then
		if [[ $download =~ "yes" ]]; then
			if [ -f "../token_header.txt" ]; then
				echo -e "OK\n"
			else
				echo -e "'token_header.txt' is missing in 'homedir'. Exiting...\n"
				rm -d ../workdir00_dataprep/ 2>/dev/null
				exit 3
			fi
		else
			if [ 0 -lt $(ls ../*fastq.gz 2>/dev/null | wc -w) ]; then
				echo -e "OK\n"
			else
				echo -e "No *.fastq.gz files in 'homedir'. Exiting...\n"
				rm -d ../workdir00_dataprep/ 2>/dev/null
				exit 3
			fi
		fi
	else
		echo -e "'renamelist.txt' is missing in 'homedir'. Exiting...\n"
		rm -d ../workdir00_dataprep/ 2>/dev/null
		exit 3
	fi
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir00_dataprep)" ]; then
		echo -e "Directory 'workdir00_dataprep' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir00_dataprep/ 2>/dev/null
		exit 3
	fi
fi
#Make folder for your data
mkdir -p $path/
mkdir 10rawreads
cd 10rawreads

#Download files from BaseSpace
if [[ $download =~ "yes" ]]; then
	echo -e "Downloading FASTQ files from BaseSpace started...\n"
	if [[ $location == "1" ]]; then
		cp /storage/$server/home/$LOGNAME/token_header.txt .
	else
		cp ../../token_header.txt .
	fi
	for (( i=$first; i<=$last; i++ )); do
		echo $i >> downloadlist.txt
	done
	if [[ $location == "1" ]]; then
		cat downloadlist.txt | parallel 'curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/{}/content -O'
	elif [[ $location == "2" ]]; then
		cat downloadlist.txt | parallel --max-procs $NSLOTS 'curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/{}/content -O'
	elif [[ $location == "0" ]]; then
		cat downloadlist.txt | parallel 'curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/{}/content -O'
	fi
	rm token_header.txt downloadlist.txt
	echo -e "\nDownloading FASTQ files from BaseSpace finished...\n"
else
	echo -e "Copying FASTQ files from home...\n"
	#Copy all *fastq.gz files from 'homedir'
	cp ../../*fastq.gz .
fi

#Copy list for file renaming & make folder structure for samples
if [[ $location == "1" ]]; then
	cp /storage/$server/home/$LOGNAME/renamelist.txt .
else
	cp ../../renamelist.txt .
fi

#Add LF at the end of last line in renamelist.txt if missing
sed -i.bak '$a\' renamelist.txt
#Delete empty lines from renamelist.txt (if any)
sed -i.bak2 '/^$/d' renamelist.txt
#Remove *.bak
rm renamelist.txt.bak renamelist.txt.bak2

for i in $(cat renamelist.txt | cut -f1)
do
	mkdir $i
done

#Rename files and move to folders according to samples (R1 and R2 files from the same sample to the same folder)
echo -e "Renaming and moving FASTQ files...\n"
cat renamelist.txt | while read -r a b
do
	mv ${b}*R1* ${a}/${a}_L001_R1_001.fastq.gz
	mv ${b}*R2* ${a}/${a}_L001_R2_001.fastq.gz
done

# Prepare samples list
cat renamelist.txt | cut -f1 > SamplesFileNames.txt
rm renamelist.txt 
cd ..
cp -r 10rawreads $path

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir00_dataprep
fi

echo -e "Script HybPhyloMaker0a finished...\n"
