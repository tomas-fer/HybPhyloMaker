#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker1a2_rawprocessParallel_summary
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 2
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker1a2_rawprocessParallel_summary
#$ -o HybPhyloMaker1a2_rawprocessParallel_summary.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *             Script 01a2 - Raw data processing in parallel summary            *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#This makes a summary table of raw data filtering after HybPhyloMaker1a_rawprocessParallel.sh

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker1a2 is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker1a2 is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir01a2
	cd workdir01a2
else
	echo -e "\nHybPhyloMaker1a2 is running locally..."
	echo -e "This summary is only for cluster environment. Exiting...\n"
	exit 3
fi

#Write log
logname=HPM1a2
echo -e "HybPhyloMaker1a2: summary of data filtering in parallel" > ${logname}.log
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
for set in data adapterfile; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

# echo -e "\nGenerating file for download and for import to Geneious...\n"
# #Copy all -all-no-dups.fas from all subfolders to folder 'for_Geneious' (for easy import to Geneious)
# mkdir $path/20filtered/for_Geneious
# find $path/20filtered/ -name '*no-dups.fastq.gz' -exec cp -t $path/20filtered/for_Geneious/ {} +
# #Pack and combine all *all-no-dups.fas files to one archive for easy download
# tar cfz $path/20filtered/for_Geneious/$data-no-dups.tar.gz -C $path/20filtered/for_Geneious/ . 2>/dev/null
# #Delete copies of fastq files in 'for Geneious' folder
# rm $path/20filtered/for_Geneious/*.fastq.gz

#Copy file list and summary logs
echo -e "Copying logs...\n"
cp $path/10rawreads/SamplesFileNames.txt .
cat SamplesFileNames.txt | while read line; do
	mkdir $line
	cp $path/20filtered/$line/*.log $line
done

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
	number=$(cut -d'_' -f2 <<< $line)
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
	withoutdups=`cat duplicate-removal-fastuniq.log | tail -n 1`
	#withoutdups=`cat duplicate-removal.log | head -n 2 | tail -n 1 | cut -d' ' -f2`
	#Calculate percentage of duplicate reads removed
	percdups=$(echo -e "scale=4;100 - 100 * ($withoutdups / $aftertrimming)" | bc)
	
	cd ..
	echo -e "$number\t$genus\t$species\t$npairs\t$nreads\t$withoutPhiX\t$percPhiX\t$bs\t$fos\t$ros\t$aftertrimming\t$perctrimm\t$withoutdups\t$percdups" >> reads_summary.txt
	rm -r $line
done

#Copy table to storage
cp reads_summary.txt $path/20filtered/

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/20filtered/

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir01a2
fi

echo -e "\nScript HybPhyloMaker1a2 finished...\n"
