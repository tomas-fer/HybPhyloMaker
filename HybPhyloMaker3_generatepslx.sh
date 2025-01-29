#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker3_generate_pslx
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker3_generate_pslx
#$ -o HybPhyloMaker3_generate_pslx.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *         Script 03 - Process consensus after mapping, make pslx files         *
# *                                   v.1.8.0c                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

# Input: Consensus sequences from HybPhyloMaker2 or Geneious: must be named consensus.fasta or consensus_cpDNA.fasta
# it is multiple fasta file with names (in case of Geneious mapping):
# e.g., Camptandra-latifolia_S4-all-no-dups_assembled_to_Curcuma_exons_reference_400Ns__consensus_sequence)
# prepared in /storage/$server/home/$LOGNAME/data/30consensus/

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker3 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add blat-suite-34
	module add r/4.4.0-gcc-10.2.1-ssuwpvb
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
	module load tools/R/3.2.1
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

#Write log
logname=HPM3
echo -e "HybPhyloMaker3: generate PSLX from consensus" > ${logname}.log
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
for set in data cp probes cpDNACDS conscall otherpslx othersource minident; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

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

#Copy fasta from home folder to scratch/workdir
cp -r $path/$type/30consensus/* .
#Make a new folder for results
mkdir -p $path/$type
mkdir $path/$type/40contigs
#Copy necessary script
cp $source/histogram.r .

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
echo -e "Generating pslx files using BLAT..."
#Copy other transcriptome/genome data from home to scratch/workdir (must be named with suffix *.fas)
if [ "$othersource" != "" ] && [ "$othersource" != "NO" ]; then 
	cp -r $othersourcepath/* .
fi

#Make a list of all files with contigs
ls *.fas > contig_names.txt

if [[ $cp =~ "yes" ]]; then
	#Copy cpDNA reference
	cp -r $source/$cpDNACDS .
else
	#Copy reference
	cp -r $source/$probes .
fi
#Make a new folder for results
mkdir $path/$type/50pslx

#A loop to process all contig files specified in contig_names.txt
for contigfile in $(cat contig_names.txt)
do
	echo -e "\nProcessing $contigfile..."
	if [[ $cp =~ "yes" ]]; then
		blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $cpDNACDS $contigfile ${contigfile}.pslx
	else
		blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $probes $contigfile ${contigfile}.pslx
		#Modify PSLX if the consensus was called using 'ConsensusFixer'
		#BLAT change all ambiguous bases to 'n', following lines put ambiguities back
		if [[ $conscall =~ "consensusfixer" ]]; then
			echo -e "Modifying PSLX..."
			head -5 ${contigfile}.pslx > pslx_header.txt #get first 5 lines as a header
			sed -i '1,5d' ${contigfile}.pslx #delete header (first 5 linees)
			cut -f22 < ${contigfile}.pslx > pslx_sequences.txt #extract column 22, i.e. column with samples-specific sequences
			sed -i "s/,$/\n/" pslx_sequences.txt #replace ',' at thee end of line by EOL (separate lines by empty line)
			sed -i "s/,/\n/g" pslx_sequences.txt #replace each ',' by EOL (put each fragment to a separate line)
			sed -i 's/n/\[ACGTRYSWKMBDHVN\]/g' pslx_sequences.txt #replace any 'n' by every possible DNA letter incl. ambiguities (as regex)
			cat pslx_sequences.txt | while read line; do
				if [ ! -z "$line" ]; then
					grep -io -m 1 $line $contigfile | head -1 >> pslx_sequences_ambig.txt
				else
					echo >> pslx_sequences_ambig.txt
				fi
			done
			#convert output to lowercase
			tr A-Z a-z < pslx_sequences_ambig.txt > pslx_sequences_ambigLower.txt
			#format back to conform PSLX standards
			sed -i '/^$/!s/$/,/' pslx_sequences_ambigLower.txt #add ',' on non-empty lines
			tr '\n' 'Q' < pslx_sequences_ambigLower.txt | sed "s/,QQ/,\n/g" | sed 's/Q//g' > pslx_sequences_ambigLowerModif.txt
			#combine with original PSLX
			cat ${contigfile}.pslx | cut -f 1-21 > first.txt #extract first 21 columns
			cat ${contigfile}.pslx | cut -f 23 > second.txt #extract 23rd column
			paste first.txt pslx_sequences_ambigLowerModif.txt second.txt > final.txt #add the modified 22nd column
			cat pslx_header.txt final.txt > ${contigfile}.pslx #add original header
			rm pslx_header.txt pslx_sequences.txt pslx_sequences_ambig.txt pslx_sequences_ambigLower.txt pslx_sequences_ambigLowerModif.txt first.txt second.txt final.txt
		fi
	fi
	cp $contigfile.pslx $path/$type/50pslx
done

#Calculate similarity to reference
echo -e "\nCalculating similarity to reference..."
#Make a list of pslx files (shorten name to genus-species_code)
ls *.pslx | cut -d'_' -f1,2 | cut -d'.' -f1 > ListOfpslx.txt
#Loop over pslx files
for pslx in $(cat ListOfpslx.txt)
do
	echo $pslx
	#Delete first five lines
	sed -i '1,5d' ${pslx}*.pslx
	#Print 1st column / ( 1st column + 2nd column ). 1st column = number of matches, 2nd column = number of mismatches
	awk '{ print ( $1 / ($1 + $2) ) * 100}' ${pslx}*.pslx > similarities
	#Call R script to make a histogram (produce PNG file called 'simil.png')
	R --slave -f histogram.r $pslx > /dev/null
	#Calculate average of similarity values
	awk '{ sum += $1 } END { print sum / NR }' similarities >> summary.txt
	cat similarities > ${pslx}_similaritytoreference.txt
	mv simil.png ${pslx}_similaritytoreference.png
done

paste ListOfpslx.txt summary.txt > similarities_summary.txt

#Make a dir for results
mkdir $path/$type/50pslx/pslxsimil
mkdir $path/$type/50pslx/pslxsimil/similarities
mkdir $path/$type/50pslx/pslxsimil/histograms
#Copy data to home
cp similarities_summary.txt $path/$type/50pslx/pslxsimil
cp *_similaritytoreference.txt $path/$type/50pslx/pslxsimil/similarities
cp *_similaritytoreference.png $path/$type/50pslx/pslxsimil/histograms

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/50pslx

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
