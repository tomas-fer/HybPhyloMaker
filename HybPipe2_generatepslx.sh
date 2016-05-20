#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=4d
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=4gb
#PBS -l scratch=8gb
#PBS -N HybPipe2_generate_pslx
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe2_generate_pslx
#$ -o HybPipe2_generate_pslx.log


# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *         Script 02 - Process consensus after mapping, make pslx files         *
# *                                   v.1.0.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#Input: Consensus sequences from Geneious: must be named consensus.fasta (it is multiple fasta file with names,
# e.g., Camptandra-latifolia_S4-all-no-dups_assembled_to_Curcuma_exons_reference_400Ns__consensus_sequence)
# prepared in /storage/$server/home/$LOGNAME/data/30consensus/

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo "Metacentrum..."
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
elif [[ $HOSTNAME == *local* ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir workdir02
	cd workdir02
	#Add necessary modules
	module load bioinformatics/blat/36x1
else
	echo "Local..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir workdir02
	cd workdir02
fi

#Copy fasta from home folder to scratch/workdir
cp -r $path/30consensus/* .
#Make a new folder for results
mkdir $path/40contigs

#-----------------------GENEIOUS CONSENSUS SEQUENCE MODIFICATION-----------------------
echo -e "\nParsing Geneious consensus sequence output...\n"
#Modify Windows EOLs to Unix EOLs (i.e., LF only)
sed -i 's/\x0D$//' consensus.fasta
#Remove trailing '?'s, remove unwanted part of the file name (everything after the second '-'), add '_consensus_sequence' to all headers
#and split multiple fasta from Geneious into individual fasta sequences
cat consensus.fasta | sed 's/^?*//' | cut -d"-" -f1,2 | sed '/>/s/.*/&_consensus_sequence/' | awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'
#Former version
#cat consensus.fasta | sed 's/^?*//' | sed "s/-all-no-dups_assembled_to_${pseudoref}_//" | awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}'
rm consensus.fasta
#Make a file with names of a all fasta files in a folder. Name is modified to consist only of first two parts of the name separated by '_' 
ls *.fasta | cut -d'_' -f 1,2 > listOfConsensusFiles.txt
#A loop for preparing assemblies (consensus sequences) of individual exons from Geneious consensus
for file in $(cat listOfConsensusFiles.txt)
do
	#Modify Windows EOLs to Unix EOLs (i.e., LF only)
	#sed -i 's/.$//' $file\_consensus_sequence.fasta
	
	#Delete first line in a fasta file, i.e. header
	sed -i 1d $file\_consensus_sequence.fasta
	
	#Replace sequence of '?' by new line (\n), put bash variable ($file) to 'val' which is available to print with awk,
	#print '>Contig'+number(NR;increased by one each step)+species name (val), then to next line print sequence
	cat $file\_consensus_sequence.fasta | tr -s '?' '\n' | awk -v val=$file '{ print ">Contig" NR "_" val "\n" $1 m}' > $file\_contigs.fas
	#Copy data from scratch to home folder
	cp $file\_contigs.fas $path/40contigs
done
echo -e "\nFinished Geneious consensus sequence modifications...\n"

#-----------------------BLAT ASSEMBLIES TO REFERENCE-----------------------
echo -e "Generating pslx files using BLAT...\n"
#Copy other transcriptome/genome data from home to scratch/workdir (must be named with suffix *.fas)
cp -r $othersourcepath/* .
#Copy reference
cp -r $source/$probes .
#Make a new folder for results
mkdir $path/50pslx

#Make a list of all files with contigs
ls *.fas > contig_names.txt

#A loop to process all contig files specified in contig_names.txt
for contigfile in $(cat contig_names.txt)
do
	blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $probes $contigfile $contigfile.pslx
	cp $contigfile.pslx $path/50pslx
done
echo -e "\nFinished BLAT...\n"

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir02
fi

echo -e "\nScript HybPipe2 finished..."
