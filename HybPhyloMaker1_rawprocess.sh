#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2d
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=4gb
#PBS -l scratch=80gb
#PBS -N HybPhyloMaker1_rawprocess
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 12
#$ -q sThC.q
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker1_rawprocess
#$ -o HybPhyloMaker1_rawprocess.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                        Script 01 - Raw data processing                       *
# *                                   v.1.2.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker1 is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add bowtie2-2.2.4
	module add samtools-0.1.19
	module add bam2fastq-1.1.0
	module add trimmomatic-0.32
	module add fastx-0.0.13
	module add jdk-7
	module add perl-5.10.1
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker1 is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir01
	cd workdir01
	#Add necessary modules
	module load bioinformatics/bam2fastq/1.1.0
	module load bioinformatics/bowtie2/2.2.6
	module load bioinformatics/fastxtoolkit/0.0.13
	module load bioinformatics/samtools/1.3
	#module load bioinformatics/trimmomatic/0.33
	module load java/1.7
else
	echo -e "\nHybPhyloMaker1 is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir01
	cd workdir01
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir01)" ]; then
		echo -e "Directory 'workdir01' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir01/ 2>/dev/null
		exit 3
	else
		if [ -d "$path/20filtered" ]; then
			echo -e "Directory '$path/20filtered' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir01/ 2>/dev/null
			exit 3
		fi
	fi
else
	if [ -d "$path/20filtered" ]; then
		echo -e "Directory '$path/20filtered' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir01/ 2>/dev/null
		exit 3
	fi
fi

#Test data structure
echo -en "Testing input data structure..."
if [ -d "$path/10rawreads" ]; then #Test if 10rawreads folder exists
	if [ -f "$path/10rawreads/SamplesFileNames.txt" ]; then #Test if SamplesFileNames.txt exists
		#Copy SamplesFileNames.txt and modify it
		cp $path/10rawreads/SamplesFileNames.txt .
		#Add LF at the end of last line in SamplesFileNames.txt if missing
		sed -i.bak '$a\' SamplesFileNames.txt
		#Delete empty lines from SamplesFileNames.txt (if any)
		sed -i.bak2 '/^$/d' SamplesFileNames.txt
		for sample in $(cat SamplesFileNames.txt); do
			if [ ! -d "$path/10rawreads/$sample" ]; then #Test if each samples-specific folder exists
				echo "Directory $sample does not exist.\n"
				rm -d ../workdir01/ 2>/dev/null
				exit 3
			else
				if [ ! -f "$path/10rawreads/$sample/${sample}_"*"R1"*".fastq.gz" ] || [ ! -f "$path/10rawreads/$sample/${sample}_"*"R2"*".fastq.gz" ]; then #Test if FASTQ.gz files exist
					echo "Proper fastq.gz files missing in $sample folder...\n"
					rm -d ../workdir01/ 2>/dev/null
					exit 3
				fi
			fi
		done
	else
		echo "List of samples (SamplesFileNames.txt) is missing. Should be in 10rawreads...\n"
		rm -d ../workdir01/ 2>/dev/null
		exit 3
	fi
else
	echo "Folder 10rawreads does not exist within your homedir.\n"
	rm -d ../workdir01/ 2>/dev/null
	exit 3
fi
echo -e "OK for running HybPhyloMaker...\n"

#Copy raw reads folders and sample list to scratch/workdir
cp -r $path/10rawreads/* .
#Copy necessary files
cp $source/fastq2fasta.pl .
cp $source/PhiX.fsa .
cp $source/NEBNext-PE.fa .
cp $source/trimmomatic-0.33.jar .
#Make a new folder for results
mkdir $path/20filtered
#Create a reference PhiX index
echo -en "Creating PhiX reference..."
bowtie2-build PhiX.fsa phiX.index 1>buildPhiX.log

#Add LF at the end of last line in SamplesFileNames.txt if missing
sed -i.bak '$a\' SamplesFileNames.txt
#Delete empty lines from SamplesFileNames.txt (if any)
sed -i.bak2 '/^$/d' SamplesFileNames.txt

#A loop to process all samples in folders named as specified in SaplesFileNames.txt
numberfiles=$(cat SamplesFileNames.txt | wc -l)
calculating=0
for file in $(cat SamplesFileNames.txt)
do
	calculating=$((calculating + 1))
	echo -e "\nProcessing sample $file ($calculating out of $numberfiles)"
	#Make a new folder for results
	mkdir $path/20filtered/$file
	#Go to sample folder
	cd $file
	#Unzip all fastq.gz files (and store only unziped files)
	echo "Unzipping..."
	gunzip *.fastq.gz
	#Get variables (read1 and read2) from list of fastq files in the folder
	read1=$(ls *.fastq | sed -n 1p)
	read2=$(ls *.fastq | sed -n 2p)
	#Map the Hyb-Seq reads of each accession to the phiX.index
	echo "Removing PhiX..."
	bowtie2 -x ../phiX.index -1 $read1 -2 $read2 -S $file.sam >bowtie2.log 2>&1
	rm  *.fastq
	#Convert .sam (resulting file of Bowtie 2) to .bam with SAMtools
	samtools view -bT ../PhiX.fsa $file.sam > $file.bam 2>/dev/null
	rm $file.sam
	#Remove the PhiX reads
	bam2fastq -o $file-noPhiX#.fq --no-aligned $file.bam >PhiX-removal.log 2>&1
	rm $file.bam
	#The resulting files need to be compressed for the next step of quality trimming
	gzip $file-noPhiX_1.fq
	gzip $file-noPhiX_2.fq
	#Adapter and general quality trimming
	echo "Adapter removal & quality trimming..."
	if [[ $location == "1" ]]; then
		java -jar /software/trimmomatic-0.32/dist/jar/trimmomatic-0.32.jar PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../NEBNext-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1
		#trimmomatic PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../NEBNext-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1
	elif [[ $location == "2" ]]; then
		java -d64 -server -XX:MaxHeapSize=10g -jar ../trimmomatic-0.33.jar PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../NEBNext-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1
	else
		java -jar ../trimmomatic-0.33.jar PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../NEBNext-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1
	fi	
	#Convert .fastq files to .fasta
	perl ../fastq2fasta.pl -a $file-1P
	perl ../fastq2fasta.pl -a $file-1U
	perl ../fastq2fasta.pl -a $file-2P
	perl ../fastq2fasta.pl -a $file-2U
	rm $file-noPhiX_1.fq.gz
	rm $file-noPhiX_2.fq.gz
	#Combine the four Trimmomatic output files in one file
	cat $file-1P.fa $file-1U.fa $file-2P.fa $file-2U.fa > $file-all.fa
	rm $file-1P* $file-1U* $file-2P* $file-2U*
	#Remove duplicate reads
	echo "Removing duplicates..."
	fastx_collapser -v -i $file-all.fa -o $file-all-no-dups.fas >duplicate-removal.log 2>&1
	#Copy results from scratch to storage
	if [[ $location == "1" ]]; then
		cp $file-all-no-dups.fas ${path}/20filtered/${file}
		cp $file-all.fa ${path}/20filtered/${file}
		cp *.log $path/20filtered/$file
	else
		cp $file-all-no-dups.fas ../${path}/20filtered/${file}
		cp $file-all.fa ../$path/20filtered/$file
		cp *.log ../$path/20filtered/$file
	fi
	rm $file-all-no-dups.fas $file-all.fa
	cd ..
done

echo -e "\nGenerating file for download and for import to Geneious...\n"
#Copy all -all-no-dups.fas from all subfolders to folder 'for_Geneious' (for easy import to Geneious)
mkdir $path/20filtered/for_Geneious
find $path/20filtered/ -name '*all-no-dups*' -exec cp -t $path/20filtered/for_Geneious/ {} +
#Pack and combine all *all-no-dups.fas files to one archive for easy download

tar cfz $path/20filtered/for_Geneious/$data-all-no-dups.tar.gz -C $path/20filtered/for_Geneious/ . 2>/dev/null
#Delete copies of fasta files in 'for Geneious' folder
rm $path/20filtered/for_Geneious/*.fas

#Summary of basic reads processing (nr. reads, PhiX filtering, quality trimming, duplicate removal)
#Produce tab-separated table
echo -e "Creating summary table..."
#Write headers (number of pairs, number of reads, number of reads after PhiX removal, percentage od PhiX reads,
#quality trimming - both surviving, forward only surviving, reverse only surviving, reads after quality trimming,
#percentage of trimmed reads, reads after duplicate removal, percentage of duplicated reads)
echo -e "Sample no.\tGenus\tSpecies\tNr. of pairs\tNr. of reads\tNr. of reads without PhiX\t% PhiX reads\tBoth surviving\tForward only surviving\tReverse only surviving\tNr. reads after quality trimming\t% quality trimmed reads\tNr. reads without duplicates\t% duplicates" > reads_summary.txt

cat SamplesFileNames.txt | while read line
do
	echo "Processing $line"
	number=$(cut -d'_' -f2 <<< $line | sed 's/S//')
	genus=$(cut -d'-' -f1 <<< $line)
	species=$(cut -d'_' -f1 <<< $line | cut -d'-' -f2)
	cd $line
	#Number of pairs: first number on first line in bowtie2.log
	npairs=`cat bowtie2.log | head -n 1 | cut -d' ' -f1`
	#Number of reads: first number on third line in PhiX-removal.log
	nreads=`cat PhiX-removal.log | head -n 3 | tail -n 1 | cut -d' ' -f1`
	#Number of reads after removal of PhiX reads: first number on fourth line in PhiX-removal.log
	withoutPhiX=`cat PhiX-removal.log | head -n 4 | tail -n 1 | cut -d' ' -f1`
	#Percentage of PhiX reads: first number on last line in bowtie2.log
	percPhiX=`cat bowtie2.log | tail -n 1 | cut -d'%' -f1`
	#Number of both surviving reads after Trimmomatic
	bs=`cat Trimmomatic.log | grep 'Input Read Pairs:' | cut -d' ' -f7`
	#Number of forward only surviving reads after Trimmomatic
	fos=`cat Trimmomatic.log | grep 'Input Read Pairs:' | cut -d' ' -f12`
	#Number of reverse only surviving reads after Trimmomatic
	ros=`cat Trimmomatic.log | grep 'Input Read Pairs:' | cut -d' ' -f17`
	# Counting number of reads after quality trimming
	let aftertrimming="(2 * $bs) + $fos + $ros"
	#Calculate percentage of quality trimmed reads
	perctrimm=$(echo -e "scale=4;100 - 100 * ($aftertrimming / $withoutPhiX)" | bc)
	#Number of reads after duplicate removal
	withoutdups=`cat duplicate-removal.log | head -n 2 | tail -n 1 | cut -d' ' -f2`
	#Calculate percentage of duplicate reads removed
	percdups=$(echo -e "scale=4;100 - 100 * ($withoutdups / $aftertrimming)" | bc)
	
	cd ..
	echo -e "$number\t$genus\t$species\t$npairs\t$nreads\t$withoutPhiX\t$percPhiX\t$bs\t$fos\t$ros\t$aftertrimming\t$perctrimm\t$withoutdups\t$percdups" >> reads_summary.txt
	rm -r $line
done

#Copy table to storage
cp reads_summary.txt $path/20filtered/

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir01
fi

echo -e "\nScript HybPhyloMaker1 finished...\n"
