#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker1a_rawprocessParallel
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker1a_rawprocessParallel
#$ -o HybPhyloMaker1a_rawprocessParallel.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                  Script 01a - Raw data processing in parallel                *
# *                                   v.1.8.0b                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#This only works on a cluster where multiple jobs are submitted

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker1a is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker1a is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir01a
	cd workdir01a
else
	echo -e "\nHybPhyloMaker1a is running locally...\n"
	echo -e "...parallel processing is not supported, run HybPhyloMaker1_rawprocess.sh instead"
	echo -e "Exiting..."
	exit 3
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir01)" ]; then
		echo -e "Directory 'workdir01' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir01a/ 2>/dev/null
		exit 3
	else
		if [ -d "$path/20filtered" ]; then
			echo -e "Directory '$path/20filtered' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir01a/ 2>/dev/null
			exit 3
		fi
	fi
else
	if [ -d "$path/20filtered" ]; then
		echo -e "Directory '$path/20filtered' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir01a/ 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM1a
echo -e "HybPhyloMaker1a: raw reads filtering in parallel" > ${logname}.log
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
for set in data adapterfile; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

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
				echo -e "Directory $sample does not exist.\n"
				rm -d ../workdir01a/ 2>/dev/null
				exit 3
			else
				if [ ! -f "$path/10rawreads/$sample/${sample}_"*"R1"*".fastq.gz" ] || [ ! -f "$path/10rawreads/$sample/${sample}_"*"R2"*".fastq.gz" ]; then #Test if FASTQ.gz files exist
					echo -e "Proper fastq.gz files missing in $sample folder...\n"
					rm -d ../workdir01a/ 2>/dev/null
					exit 3
				fi
			fi
		done
	else
		echo "List of samples (SamplesFileNames.txt) is missing. Should be in 10rawreads...\n"
		rm -d ../workdir01a/ 2>/dev/null
		exit 3
	fi
else
	echo -e "Folder 10rawreads does not exist within your homedir.\n"
	rm -d ../workdir01a/ 2>/dev/null
	exit 3
fi
echo -e "OK for running HybPhyloMaker...\n"

#Make a new folder for results
mkdir -p $path/20filtered
#Create new 'submitRawProcessJobs.sh' and make it executable
touch $path/20filtered/submitRawProcessJobs.sh
chmod +x $path/20filtered/submitRawProcessJobs.sh

echo -en "\nGenerating multiple jobs (1 job per sample)..."

#Add LF at the end of last line in SamplesFileNames.txt if missing
sed -i.bak '$a\' SamplesFileNames.txt
#Delete empty lines from SamplesFileNames.txt (if any)
sed -i.bak2 '/^$/d' SamplesFileNames.txt

for file in $(cat SamplesFileNames.txt)
do
	echo '#!/bin/bash' >> ${file}.sh
	echo '#----------------MetaCentrum----------------' >> ${file}.sh
	echo '#PBS -l walltime=2:0:0' >> ${file}.sh
	echo '#PBS -l select=1:ncpus=2:mem=16gb:scratch_local=24gb' >> ${file}.sh
	echo '#PBS -j oe' >> ${file}.sh
	echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> ${file}.sh
	echo '#PBS -N rawprocess_for_'"${file}" >> ${file}.sh
	echo '#-------------------HYDRA-------------------' >> ${file}.sh
	echo '#$ -S /bin/bash' >> ${file}.sh
	echo '#$ -pe mthread 2' >> ${file}.sh
	echo '#$ -q mThC.q' >> ${file}.sh
	echo '#$ -l mres=4G,h_data=4G,h_vmem=4G' >> ${file}.sh
	echo '#$ -cwd' >> ${file}.sh
	echo '#$ -j y' >> ${file}.sh
	echo '#$ -N rawprocess_for_'"${file}" >> ${file}.sh
	echo '#$ -o rawprocess_for_'"${file}"'.log' >> ${file}.sh
	echo '#Complete path and set configuration for selected location' >> ${file}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${file}.sh
	echo '  #Add necessary modules' >> ${file}.sh
	echo '  module add bowtie2-2.2.4' >> ${file}.sh
	echo '  module add samtools-0.1.19' >> ${file}.sh
	echo '  module add bam2fastq-1.1.0' >> ${file}.sh
	echo '  module add trimmomatic-0.32' >> ${file}.sh
	echo '  module add fastuniq-1.1' >> ${file}.sh
	echo '  module add jdk-7' >> ${file}.sh
	echo '  module add perl-5.10.1' >> ${file}.sh
	echo '  cd $SCRATCHDIR' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  #Add necessary modules' >> ${file}.sh
	echo '  module load bioinformatics/bowtie2/2.2.9' >> ${file}.sh
	echo '  module load bioinformatics/samtools/1.3' >> ${file}.sh
	echo '  module load bioinformatics/bam2fastq/1.1.0' >> ${file}.sh
	echo '  module load bioinformatics/fastuniq/1.1' >> ${file}.sh
	echo '  module load java/1.7' >> ${file}.sh
	echo '  mkdir workdir01a_'"${file}" >> ${file}.sh
	echo '  cd workdir01a_'"${file}" >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo 'path='"$path" >> ${file}.sh
	echo 'source='"$source" >> ${file}.sh
	echo 'file='"$file" >> ${file}.sh
	echo 'location='"$location" >> ${file}.sh
	echo 'adapterfile='"$adapterfile" >> ${file}.sh
	echo '#Copy raw reads folders and sample list to scratch/workdir' >> ${file}.sh
	echo 'cp -r $path/10rawreads/'"${file}"' .' >> ${file}.sh
	echo '#Copy necessary files' >> ${file}.sh
	echo 'cp $source/fastq2fasta.pl .' >> ${file}.sh
	echo 'cp $source/PhiX.fsa .' >> ${file}.sh
	echo 'cp $source/'"${adapterfile}"' .' >> ${file}.sh
	echo 'cp $source/trimmomatic-0.33.jar .' >> ${file}.sh
	echo '#Create a reference PhiX index' >> ${file}.sh
	echo 'echo -en "Creating PhiX reference..."' >> ${file}.sh
	echo 'bowtie2-build PhiX.fsa phiX.index 1>buildPhiX.log' >> ${file}.sh
	echo '#Make a new folder for results' >> ${file}.sh
	echo 'mkdir $path/20filtered/$file' >> ${file}.sh
	echo '#Go to sample folder' >> ${file}.sh
	echo 'cd $file' >> ${file}.sh
	echo '#Unzip all fastq.gz files (and store only unziped files)' >> ${file}.sh
	echo 'echo "Unzipping..."' >> ${file}.sh
	echo 'gunzip *.fastq.gz' >> ${file}.sh
	echo '#Get variables (read1 and read2) from list of fastq files in the folder' >> ${file}.sh
	echo 'read1=$(ls *.fastq | sed -n 1p)' >> ${file}.sh
	echo 'read2=$(ls *.fastq | sed -n 2p)' >> ${file}.sh
	echo '#Map the Hyb-Seq reads of each accession to the phiX.index' >> ${file}.sh
	echo 'echo "Removing PhiX..."' >> ${file}.sh
	echo 'bowtie2 -x ../phiX.index -1 $read1 -2 $read2 -S $file.sam >bowtie2.log 2>&1' >> ${file}.sh
	echo 'rm  *.fastq' >> ${file}.sh
	echo '#Convert .sam (resulting file of Bowtie 2) to .bam with SAMtools' >> ${file}.sh
	echo 'samtools view -bT ../PhiX.fsa $file.sam > $file.bam 2>/dev/null' >> ${file}.sh
	echo 'rm $file.sam' >> ${file}.sh
	echo '#Remove the PhiX reads' >> ${file}.sh
	echo 'bam2fastq -o $file-noPhiX#.fq --no-aligned $file.bam >PhiX-removal.log 2>&1' >> ${file}.sh
	echo 'rm $file.bam' >> ${file}.sh
	echo '#The resulting files need to be compressed for the next step of quality trimming' >> ${file}.sh
	echo 'gzip $file-noPhiX_1.fq' >> ${file}.sh
	echo 'gzip $file-noPhiX_2.fq' >> ${file}.sh
	echo '#Adapter and general quality trimming' >> ${file}.sh
	echo 'echo "Adapter removal & quality trimming..."' >> ${file}.sh
	echo 'if [[ $location == "1" ]]; then' >> ${file}.sh
	echo '  java -jar /software/trimmomatic-0.32/dist/jar/trimmomatic-0.32.jar PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../${adapterfile}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1' >> ${file}.sh
	echo '  #trimmomatic PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../NEBNext-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1' >> ${file}.sh
	echo 'elif [[ $location == "2" ]]; then' >> ${file}.sh
	echo '  java -d64 -server -XX:MaxHeapSize=4g -jar ../trimmomatic-0.33.jar PE -phred33 $file-noPhiX_1.fq.gz $file-noPhiX_2.fq.gz $file-1P $file-1U $file-2P $file-2U ILLUMINACLIP:../${adapterfile}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36 >Trimmomatic.log 2>&1' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo 'rm $file-noPhiX_1.fq.gz' >> ${file}.sh
	echo 'rm $file-noPhiX_2.fq.gz' >> ${file}.sh
	echo '#Remove duplicate reads' >> ${file}.sh
	echo 'echo "Removing duplicates..."' >> ${file}.sh
	echo 'ls $file-?P > fastqlist.txt' >> ${file}.sh
	echo 'fastuniq -i fastqlist.txt -t q -o ${file}-1P_no-dups.fastq -p ${file}-2P_no-dups.fastq' >> ${file}.sh
	echo '#Count number of duplicates' >> ${file}.sh
	echo '#number of unpaired reads' >> ${file}.sh
	echo 'u1=`expr $(cat ${file}-1U | wc -l) / 4`' >> ${file}.sh
	echo 'u2=`expr $(cat ${file}-2U | wc -l) / 4`' >> ${file}.sh
	echo '#paired-read number before duplicate removal' >> ${file}.sh
	echo 'a=`expr $(cat ${file}-1P | wc -l) / 4`' >> ${file}.sh
	echo 'echo "Nr. of reads before duplicate removal:" >> duplicate-removal-fastuniq.log' >> ${file}.sh
	echo 'before=`expr $a + $a + $u1 + $u2`' >> ${file}.sh
	echo 'echo $before >> duplicate-removal-fastuniq.log' >> ${file}.sh
	echo '#paired-read number after duplicate removal' >> ${file}.sh
	echo 'echo "Nr. of reads after duplicate removal:" >> duplicate-removal-fastuniq.log' >> ${file}.sh
	echo 'b=`expr $(cat ${file}-1P_no-dups.fastq | wc -l) / 4`' >> ${file}.sh
	echo 'after=`expr $b + $b + $u1 + $u2`' >> ${file}.sh
	echo 'echo $after >> duplicate-removal-fastuniq.log' >> ${file}.sh
	echo '#gzip trimmed and deduplicated reads' >> ${file}.sh
	echo 'gzip ${file}-1P' >> ${file}.sh
	echo 'gzip ${file}-1U' >> ${file}.sh
	echo 'gzip ${file}-2P' >> ${file}.sh
	echo 'gzip ${file}-2U' >> ${file}.sh
	echo 'gzip ${file}-1P_no-dups.fastq' >> ${file}.sh
	echo 'gzip ${file}-2P_no-dups.fastq' >> ${file}.sh
	echo '#Copy results from scratch to storage' >> ${file}.sh
	echo 'if [[ $location == "1" ]]; then' >> ${file}.sh
	echo '  cp ${file}-1P_no-dups.fastq.gz ${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-2P_no-dups.fastq.gz ${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-1P.gz ${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-1U.gz ${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-2P.gz ${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-2U.gz ${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp *.log $path/20filtered/$file' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  cp ${file}-1P_no-dups.fastq.gz ../${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-2P_no-dups.fastq.gz ../${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-1P.gz ../${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-1U.gz ../${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-2P.gz ../${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp ${file}-2U.gz ../${path}/20filtered/${file}' >> ${file}.sh
	echo '  cp *.log ../$path/20filtered/$file' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo 'mv ${file}-1U.gz ${file}-1U_no-dups.fastq' >> ${file}.sh
	echo 'mv ${file}-2U.gz ${file}-2U_no-dups.fastq' >> ${file}.sh
	echo 'cd ..' >> ${file}.sh
	echo '#Clean scratch/work directory' >> ${file}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${file}.sh
	echo '  #delete scratch' >> ${file}.sh
	echo '  rm -rf $SCRATCHDIR/*' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  cd ..' >> ${file}.sh
	echo '  rm -r workdir01a_'"${file}" >> ${file}.sh
	echo 'fi' >> ${file}.sh
	
	chmod +x ${file}.sh
	if [[ $location == "1" ]]; then
		cp ${file}.sh $path/20filtered
		#qsub ${file}.sh
		echo 'qsub '"${file}"'.sh' >> $path/20filtered/submitRawProcessJobs.sh
	else
		cp ${file}.sh $path/20filtered
		cp ${file}.sh ..
		echo 'qsub '"${group}"'.sh' >> ../submitRawProcessJobs.sh
	fi
done

echo -e "done\n"

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/20filtered

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir01a
fi

echo -e "\nHybPhyloMaker 1a finished..."
if [[ $location == "2" ]]; then
	echo -e "\nGo to homedir and run submitRawProcessJobs.sh...\n"
else
	echo -e "\nGo to $path/20filtered and run submitRawProcessJobs.sh..."
	echo -e "This starts parallel filtering of raw reads."
	echo -e "\nAfter all jobs finish run script HybPhyloMaker1a2 in order to summarize the filtered data...\n"
fi
