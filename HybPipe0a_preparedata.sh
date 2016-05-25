#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=4:ppn=6
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=8gb
#PBS -N HybPipe0a_datadownloadprepare
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 24
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe0a_datadownloadprepare
#$ -o HybPipe0a_datadownloadprepare.log

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                     Script 00a - Download & prepare data                     *
# *                                   v.1.0.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


# Download fastq files from BaseSpace (optional), rename them and move to sample folders
# i.e., prepare raw reads for running HybPipe
# Works for PE reads (i.e., 2 fastq files per sample) in the following format:
# samplename_sampleID_L001_R1_001.fastq.gz
# samplename_sampleID_L001_R2_001.fastq.gz
#
# Needs two files (token_header.txt, renamelist.txt) to be in the home folder (at desired data server)


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
#                                                                                                                   #
# 4. prepare file renamelist.txt with two columns (desired sample name and 'Sample_Name' from BaseSpace)            #
#   (Sample_Name is first part of file name at BaseSpace), e.g.                                                     #
# genus1-species1_S001	Z001                                                                                        #
# genus1-species2_S002	Z002                                                                                        #
# etc.                                                                                                              #
# To ensure smooth processing with subsequent scripts in HybPipe the file names must follow the convention          #
# genus-species_code                                                                                                #
#                                                                                                                   #
#####################################################################################################################

if [[ $PBS_O_HOST == *".cz" ]]; then
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
elif [[ $HOSTNAME == *local* ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir00_ref
	cd workdir00_ref
	#Add necessary modules
	module load tools/gnuparallel/20160422
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir00_dataprep
	cd workdir00_dataprep
fi

#Make folder for your data
mkdir -p $path/
mkdir 10rawreads
cd 10rawreads

#Download files from BaseSpace
if [[ $download =~ "yes" ]]; then
	if [[ $location == "1" ]]; then
		cp /storage/$server/home/$LOGNAME/token_header.txt .
	else
		cp ../../token_header.txt .
	fi
	for (( i=$first; i<=$last; i++ )); do
		echo $i >> downloadlist.txt
	done
	cat downloadlist.txt | parallel 'curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/{}/content -O'
	rm token_header.txt downloadlist.txt
fi

#Copy all *fastq.gz files from 'homedir'
cp ../../*fastq.gz .

#Copy list for file renaming & make folder structure for samples
if [[ $location == "1" ]]; then
	cp /storage/$server/home/$LOGNAME/renamelist.txt .
else
	cp ../../renamelist.txt .
fi

for i in $(cat renamelist.txt | cut -f1)
do
	mkdir $i
done

#Rename files and move to folders according to samples (R1 and R2 files from the same sample to the same folder)
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
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir00_dataprep
fi
