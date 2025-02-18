#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:0:0
#PBS -l select=1:ncpus=12:mem=16gb:scratch_local=400gb
#PBS -j oe
#PBS -N HybPhyloMaker0a_datadownloadprepare
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker0a_datadownloadprepare
#$ -o HybPhyloMaker0a_datadownloadprepare.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                     Script 0a - Download & prepare data                      *
# *                                   v.1.8.0d                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


# Download fastq files from SRA (Short Read Archive, NCBI), ENA (European Nucleotide Archive, EBI) or BaseSpace (optional)
# rename them and move to sample folders
# Or just rename fastq.gz samples in 'home' directory and move them to sample folders
# i.e., prepare raw reads for running HybPhyloMaker
# Controled by 'download=' option in settings.cfg (yes/no/sra)

# Works for PE reads only (i.e., 2 fastq files per sample), samples are renamed to the following format:
# samplename_sampleID_L001_R1_001.fastq.gz
# samplename_sampleID_L001_R2_001.fastq.gz

# If only renaming&moving of the samples is requested (download=no), there has to be the 'renamelist.txt' in the 'home' folder
# Samples should contain *R1* or *R2* in the name!!!

# For SRA/ENA download (download=sra) - needs file 'SRRandERRlist.txt' in HybSeqSource (simple list of SRR/ERR codes, one per line)
# Does not test validity of SRR code (will be implemented in the future?)
# No 'renamelist.txt' is required

# For BaseSpace download (download=yes)
# Needs two files (token_header.txt and renamelist.txt) to be in the 'home' folder
# (at desired data server if running on cluster)
# Download from BaseSpace is parallelized (via GNU parallel)

#####################################################################################################################
# Usage (only if you require download from BaseSpace, otherwise continue with step 4):                              #
# 1. You must obtain 'token' from Illumina BaseSpace (how to do this see steps 1-5 at                               #
# https://support.basespace.illumina.com/knowledgebase/articles/403618-python-run-downloader )                      #
# 2. Save this token to text file (token_header.txt) with one line text:                                            #
# header = "x-access-token: <your-token-here>"                                                                      #
# 3. Login to BaseSpace via web browser and get ID for the project                                                  #
# (all fastq.gz files from the project will be downloaded)                                                          #
# How to do this: go to (via clicking) Projects -> <project-name> -> Samples -> <sample-name> -> <file>.fastq.gz    #
# Look at the address which should looks like                                                                       #
# https://basespace.illumina.com/projects/2016978377                                                                #
# desired ID is the last number (put this to settings.cfg)                                                          #
# 4. prepare file renamelist.txt with two columns (desired sample name and 'Sample_Name' from BaseSpace)            #
#   (Sample_Name is first unique part of file name at BaseSpace), e.g.                                              #
# genus1-species1_S001	Z001                                                                                        #
# genus1-species2_S002	Z002                                                                                        #
# etc.                                                                                                              #
# To ensure smooth processing with subsequent scripts in HybPhyloMaker the file names must follow the convention    #
# genus-species_code                                                                                                #
# i.e., for samples containing '-' in the name and downloaded from SRA this needs to be corrected manually          #
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
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	homedir=/storage/$server/home/$LOGNAME
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add parallel
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker0a is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	homedir=..
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
	homedir=..
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir00_dataprep
	cd workdir00_dataprep
fi

#Download from SRA/ENA (based on the 'SRRandERRlist.txt' in HybSeqSource)
if [[ $download =~ "sra" ]]; then
	if [ -f "$source/SRRandERRlist.txt" ]; then
		if [[ $PBS_O_HOST == *".cz" ]]; then
			module add sratoolkit
			module add entrezdirect
			module add py-ffq
			ftp-cp ftp.ncbi.nlm.nih.gov /entrez/entrezdirect xtract.Linux.gz
			gunzip -f xtract.Linux.gz
			chmod +x xtract.Linux
			mkdir -p ~/.ncbi
			cp ${homedir}/.ncbi/user-settings.mkfg ~/.ncbi
		fi
		#Make dir for results
		mkdir -p $path
		mkdir 10rawreads
		#Loop over samples, download, rename, move to folder
		echo -e "Downloading samples from SRA...\n"
		#Download SRRandERRlist.txt
		cp ${source}/SRRandERRlist.txt .
		#Remove extrace spaces from SRRandERRlist.txt
		sed -i 's/ //g' SRRandERRlist.txt
		for xrr in $(cat SRRandERRlist.txt); do
			if grep -q "^SRR[0-9]*$" <<< $xrr; then
				#get organism name
				if [[ $PBS_O_HOST == *".cz" ]]; then
					name=$(esearch -db sra -query ${xrr} | esummary | ./xtract.Linux -pattern DocumentSummary -element Organism@ScientificName Run@acc | sed 's/ /-/' | sed 's/\t/_/')
				else
					name=$(esearch -db sra -query ${xrr} | esummary | xtract -pattern DocumentSummary -element Organism@ScientificName Run@acc | sed 's/ /-/' | sed 's/\t/_/')
				fi
				if [[ -z $name ]]; then
					echo -e "${xrr} did not return valid accession in SRA...\n"
				else
					echo -e "${name}"
					#download fastq.gz files
					echo -e "downloading"
					fastq-dump --split-3 --gzip ${xrr}
					#unzip
					echo -e "unzipping"
					gzip -d ${xrr}_{1,2}.fastq.gz
					#remove everything after '+SRR' on the beginning of line, i.e. make standard fastq format with just '+' on third line
					echo -e "modifying"
					sed -i 's/^+SRR.*/+/' ${xrr}_{1,2}.fastq
					#gzip
					echo -e "gzipping"
					gzip ${xrr}_{1,2}.fastq
					#rename files
					mv ${xrr}_1.fastq.gz ${name}_L001_R1_001.fastq.gz
					mv ${xrr}_2.fastq.gz ${name}_L001_R2_001.fastq.gz
					#create a dir and move the fastq.gz files
					mkdir 10rawreads/${name}
					mv ${name}_L001_R{1,2}_001.fastq.gz 10rawreads/${name}
					echo
				fi
			elif grep -q "^ERR[0-9]*$" <<< $xrr; then
				#get organism name
				name=$(n=$(ffq ${xrr} 2>/dev/null | grep ERS | cut -d'"' -f4) && ffq ${n} 2>/dev/null | grep org | head -n1 | cut -d'"' -f4 | sed 's/ /-/')
				if [[ -z $name ]]; then
					echo -e "${xrr} did not return valid accession in ENA...\n"
				else
					fname=${name}_${xrr}
					echo -e "${fname}"
					#download fastq.gz files
					echo -e "downloading"
					wget $(ffq --ftp ${xrr} 2>/dev/null | grep '"url"' | cut -d '"' -f4) 2>/dev/null
					#rename files
					echo -e "renaming"
					mv ${xrr}_1.fastq.gz ${fname}_L001_R1_001.fastq.gz
					mv ${xrr}_2.fastq.gz ${fname}_L001_R2_001.fastq.gz
					#create a dir and move the fastq.gz files
					mkdir 10rawreads/${fname}
					mv ${fname}_L001_R{1,2}_001.fastq.gz 10rawreads/${fname}
					echo
				fi
			else
				echo -e "${xrr} is not a valid SRR/ERR code...\n"
			fi
		done
		
		#Create list of samples
		ls 10rawreads > 10rawreads/SamplesFileNames.txt
		#Copy data home
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
		
		exit 0
	else
		echo -e "List of SRR/ERR numbers for download (SRRandERRlist.txt) is missing. Should be in HybSeqSource...\n"
		rm -d ../workdir00_dataprep/ 2>/dev/null
		exit 3
	fi
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -d "$path/10rawreads" ]; then
	echo -e "Directory '$path/10rawreads' already exists. Exiting...\n"
	rm -d ../workdir00_dataprep/ 2>/dev/null
	exit 3
else
	if [ -f "$homedir/renamelist.txt" ]; then
		if [[ $download =~ "yes" ]]; then
			if [ -f "$homedir/token_header.txt" ]; then
				echo -e "OK\n"
			else
				echo -e "'token_header.txt' is missing in 'homedir'. Exiting...\n"
				rm -d ../workdir00_dataprep/ 2>/dev/null
				exit 3
			fi
		else
			if [ 0 -lt $(ls $homedir/*fastq.gz 2>/dev/null | wc -w) ]; then
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

#Download files from BaseSpace
if [[ $download =~ "yes" ]]; then
	echo -e "Downloading FASTQ files from BaseSpace requested...\n"
	#Get token
	if [[ $location == "1" ]]; then
		cp /storage/$server/home/$LOGNAME/token_header.txt .
	else
		cp ../../token_header.txt .
	fi
	
	#Get info about samples
	echo -e "Getting info about samples in the project with ID: ${projectID}..."
	curl -L -J --config ./token_header.txt ${bsserver}/v1pre3/projects/${projectID}/samples?Limit=1000 2>/dev/null > JSONproject.txt
	#Test if there is permit to access the project
	if [ ! -z "$(grep ErrorCode JSONproject.txt)" ]; then
		echo -e "You probably does not have access rights to the project. Exiting...\n"
		rm -d ../workdir00_dataprep/ 2>/dev/null
		exit 3
	fi
	#Extract sample numbers
	grep -Po '"Status":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $4}' > statusList.txt #whether Aborted or Complete
	grep -Po '"Href":.*?[^\\]",' JSONproject.txt | grep "/samples" | awk -F\" '{print $4}'| awk -F\/ '{print $3}' > samplesList.txt
	grep -Po '"SampleId":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $4}' > sampleID.txt
	grep -Po '"LibraryName":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $4}' > libName.txt
	grep -Po '"TotalReadsPF":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > readsPF.txt
	expName=$(grep -Po '"ExperimentName":.*?[^\\]"' JSONproject.txt | awk -F\" '{print $4}' | head -n1)
	#Make samples table
	echo -e "Status\tSampleID\tName\tReadsPF\tBaseSpaceID" > sampleTable.txt
	paste statusList.txt sampleID.txt libName.txt readsPF.txt samplesList.txt | grep Complete >> sampleTable.txt #remove aborted samples
	#Make samplesList without aborted
	paste statusList.txt samplesList.txt | grep Complete | awk '{ print $2 }' > tmp
	mv tmp samplesList.txt
	rm statusList.txt sampleID.txt libName.txt readsPF.txt
	echo "There are" `cat samplesList.txt | wc -l` "samples in the project '$expName'"
	#Get file IDs
	echo -e "\nGetting info about files in the project ${projectID}..."
	for i in $(cat samplesList.txt); do
		#download information about files for particular sample
		curl -L -J --config ./token_header.txt ${bsserver}/v1pre3/samples/${i}/files?Extensions=gz 2>/dev/null > JSONsamples.txt
		#get 'Id', display only them and add it to the IDs list
		grep -Po '"Id":.*?[^\\]",' JSONsamples.txt | awk -F\" '{print $4}' >> filesList.txt
		#get sizes
		grep -Po '"Size":.*?[^\\]"' JSONsamples.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> filesSize.txt
		#get file names
		grep -Po '"Path":.*?[^\\]"' JSONsamples.txt | awk -F\" '{print $4}' >> filesNames.txt
	done
	echo "There are" `cat filesList.txt | wc -l` "files"
	#Make table
	echo -e "FileName\tSize\tBaseSpaceID" > fileTable.txt
	paste filesNames.txt filesSize.txt filesList.txt >> fileTable.txt
	rm filesNames.txt filesSize.txt
	#Download individual files using parallel
	echo -e "\nDownloading fastq.gz files...\n"
	if [[ $location == "1" ]]; then
		cat filesList.txt | parallel 'curl -L -J --config token_header.txt ${bsserver}/v1pre3/files/{}/content -O'
	elif [[ $location == "2" ]]; then
		cat filesList.txt | parallel --max-procs $NSLOTS 'curl -L -J --config token_header.txt ${bsserver}/v1pre3/files/{}/content -O'
	elif [[ $location == "0" ]]; then
		cat filesList.txt | parallel 'curl -L -J --config token_header.txt ${bsserver}/v1pre3/files/{}/content -O'
	fi
	#Check file sizes of downloaded files (if they match sizes stated by BaseSpace)
	#Download incorrectly downloaded files
	echo -e "\nChecking whether file sizes are correct..."
	cat fileTable.txt | sed '1d' | while read line; do
		fileName=$(awk '{ print $1 }' <<< $line) #file name
		fileSizeDown=$(stat -c %s `awk '{ print $1 }' <<< $line`) #file size of the downloaded file
		fileSizeBS=$(awk '{ print $2 }' <<< $line) #file size extracted from BaseSpace JSON
		fileBS=$(awk '{ print $3 }' <<< $line) #file BaseSpace ID
		if [[ ${fileSizeDown} -eq ${fileSizeBS} ]]; then
			echo ${fileName} size OK
		else
			echo -e "${fileName} size incorrect. Downloading again...\n"
			mv ${fileName} ${fileName}.bak 2>/dev/null #make a backup of the wrongly downloaded file
			until [[ $(stat -c %s `awk '{ print $1 }' <<< $line` 2>/dev/null) -eq ${fileSizeBS} ]]; do
				curl -L -J --config token_header.txt ${bsserver}/v1pre3/files/${fileBS}/content -O
				echo
			done
		fi
	done
	echo -e "\nDownloading of FASTQ files from BaseSpace finished...\n"
	#Copy results home
	mkdir -p $path/00downloadinfo
	cp JSON*.txt $path/00downloadinfo
	cp sampleTable.txt $path/00downloadinfo
	cp fileTable.txt $path/00downloadinfo
	rm filesList.txt fileTable.txt JSON*.txt samplesList.txt sampleTable.txt
	
	#Download informations about the run (only if runID is set)
	if [ ! -z $runID ]; then
		#Run details (total yield, number of clusters, total reads, total reads PF...)
		echo "Getting info about run ${runID}"
		#Get info about the run
		curl -L -J --config ./token_header.txt ${bsserver}/v1pre3/runs/${runID} 2>/dev/null > JSONrun.txt
		#Test if there is permit to access the project
		if [ ! -z "$(grep ErrorCode JSONrun.txt)" ]; then
			echo -e "You probably does not have access rights to the run. Skiping run statistics...\n"
		else
			grep -Po '"ExperimentName":.*?[^\\]"' JSONrun.txt | head -n1 | awk -F\" '{print $4}' | sed 's/[:,]//g' > rundata.txt
			grep -Po '"PlatformName":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $4}' | sed 's/[:,]//g' >> rundata.txt
			grep -Po '"YieldTotal":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
			grep -Po '"Clusters":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
			grep -Po '"ClustersPf":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,}]//g' >> rundata.txt
			grep -Po '"PercentPf":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
			grep -Po '"PercentGtQ30":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
			grep -Po '"PercentGtQ30R1":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
			grep -Po '"PercentGtQ30R2":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
			#make a table
			echo -e "ExperimentName\nPlatformName\nTotalYield[Gbp]\nClusters\nClustersPF\nPercentPF\nPercentGtQ30\nPercentGtQ30R1\nPercentGtQ30R2" > runHeader.txt
			paste runHeader.txt rundata.txt > runTable.txt
			echo -e "Summary of the run $runID is in 'runTable.txt'\n"
			cp runTable.txt $path/00downloadinfo
			cp JSONrun.txt $path/00downloadinfo
			rm rundata.txt runHeader.txt JSONrun.txt runTable.txt
		fi
	fi
else
	echo -e "Copying FASTQ files from home...\n"
	#Copy all *fastq.gz files from 'homedir'
	cat renamelist.txt | while read -r a b; do
		if [[ $location == "1" ]]; then
			cp $homedir/*${b}* .
		else
			cp ../../*${b}* .
		fi
	done
fi

#Remove token
if [[ $download =~ "yes" ]]; then
	rm token_header.txt
fi

for i in $(cat renamelist.txt | cut -f1)
do
	mkdir $i
done

#Rename files and move to folders according to samples (R1 and R2 files from the same sample to the same folder)
echo -e "Renaming and moving FASTQ files...\n"
cat renamelist.txt | while read -r a b
do
	mv *${b}*R1* ${a}/${a}_L001_R1_001.fastq.gz
	mv *${b}*R2* ${a}/${a}_L001_R2_001.fastq.gz
done

# Prepare samples list
cat renamelist.txt | cut -f1 > SamplesFileNames.txt
if [[ $download =~ "yes" ]]; then
	cp renamelist.txt $path/00downloadinfo
else
	cp renamelist.txt $path
fi
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
