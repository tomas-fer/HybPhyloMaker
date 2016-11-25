#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2d
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=4gb
#PBS -l scratch=80gb
#PBS -N HybPhyloMaker1b_readmapping
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 12
#$ -q sThC.q
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker1b_readmapping
#$ -o HybPhyloMaker1b_readmapping.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                     Script 02 - Read mapping using bowtie2                   *
# *                                   v.1.3.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker1b is running on MetaCentrum...\n"
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
	#module add bcftools-1.3.1
	module add samtools-1.3
	module add perl-5.10.1
	module add gcc-4.8.4
	module add ococo-2016-11
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker1b is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir02
	cd workdir02
	#Add necessary modules
	module load bioinformatics/bowtie2/2.2.6
	module load bioinformatics/samtools/1.3
else
	echo -e "\nHybPhyloMaker1b is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir02
	cd workdir02
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir02)" ]; then
		echo -e "Directory 'workdir02' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir01/ 2>/dev/null
		exit 3
	else
		if [ -d "$path/20filtered" ]; then
			echo -e "Directory '$path/21filtered' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir01/ 2>/dev/null
			exit 3
		fi
	fi
else
	if [ -d "$path/20filtered" ]; then
		echo -e "Directory '$path/21filtered' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir01/ 2>/dev/null
		exit 3
	fi
fi

#Copy pseudoreference
probes=$(echo $probes | cut -d"." -f1)
cp $source/${probes}_with${nrns}Ns_beginend.fas .

#Make index from pseudoreference
if [[ $mapping =~ "yes" ]]; then
	bowtie2-build ${probes}_with${nrns}Ns_beginend.fas pseudoreference.index
fi
#Make a new folder for results
mkdir $path/21mapped

#Copy list of samples
cp $path/10rawreads/SamplesFileNames.txt .

#A loop to process all samples in folders named as specified in SaplesFileNames.txt
numberfiles=$(cat SamplesFileNames.txt | wc -l)
calculating=0
for file in $(cat SamplesFileNames.txt); do
	calculating=$((calculating + 1))
	echo -e "Processing sample $file ($calculating out of $numberfiles)"
	#set parameters for mapping
	score=G,20,8
	#sensitive mapping
	if [[ $mapping =~ "yes" ]]; then
		#copy fastq files
		cp $path/20filtered/${file}/${file}-1P_no-dups.fastq .
		cp $path/20filtered/${file}/${file}-2P_no-dups.fastq .
		cp $path/20filtered/${file}/${file}-1U .
		cp $path/20filtered/${file}/${file}-2U .
		echo "mapping"
		#Bowtie2 parameters are derived from --very-sensitive-local
		bowtie2 --local -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --score-min $score -x pseudoreference.index -1 ${file}-1P_no-dups.fastq  -2 ${file}-2P_no-dups.fastq  -U ${file}-1U,${file}-2U -S ${file}.sam 2>${file}_bowtie2_out.txt
		#create SAM from BAM
		echo "converting to BAM"
		samtools view -bS -o ${file}.bam ${file}.sam
		#copy results to home
		cp ${file}.bam $path/21mapped
		cp ${file}_bowtie2_out.txt $path/21mapped
	else
		echo "copying BAM"
		cp $path/21mapped/${file}.bam .
	fi
	#CONSENSUS USING OCOCO
	echo "making consensus with OCOCO"
	ococo -i ${file}.bam -c $mincov -F ${file}.fasta
	#change name in fasta file
	sed -i '1d' ${file}.fasta #delete first line
	sed -i "1i >$file" ${file}.fasta #insert fasta header as a first line
	#Remove line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${file}.fasta > tmp && mv tmp ${file}.fasta
	#put $nrns Ns to variable 'a' and $nrns ?s to variable 'b'
	a=$(printf "%0.sN" $(seq 1 $nrns))
	b=$(printf "%0.s?" $(seq 1 $nrns))
	#replace all Ns separating exons by '?'
	sed -i "s/$a/$b/g" ${file}.fasta
	#copy results to home
	cp ${file}.fasta $path/21mapped
	#delete BAM
	rm ${file}.bam
done

#Combine all fasta file into one
cat *.fasta > consensus.fasta
#Remove line breaks from fasta file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' consensus.fasta > tmp && mv tmp consensus.fasta
mkdir $path/30consensus
cp consensus.fasta $path/30consensus

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir02
fi

echo -e "\nScript HybPhyloMaker2 finished...\n"
