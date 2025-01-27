#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker12_exclude_samples
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker12_exclude_samples
#$ -o HybPhyloMaker12_exclude_samples.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                         Script 12 - Exclude samples                          *
# *                                   v.1.8.0d                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Input: $type/60mafft (or 61mafft_corrected) and 'excludelist.txt' in HybSeqSource
#Output: $type/$selection/60mafft (or 61mafft_corrected) and $type/$selection/70concatenated_exon_alignments (or 80concatenated_exon_alignments_corrected)
#with exon and loci alignments without samples defined in 'excludelist.txt'

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker12 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	otherpslxcppath=/storage/$server/home/$LOGNAME/$otherpslxcp
	#Add necessary modules
	#module add python-2.7.6-gcc
	#module add python-2.7.6-intel
	#module add debian8-compat
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker12 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	otherpslxcppath=../$otherpslxcp
	#Make and enter work directory
	mkdir -p workdir12
	cd workdir12
else
	echo -e "\nHybPhyloMaker12 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	otherpslxcppath=../$otherpslxcp
	#Make and enter work directory
	mkdir -p workdir12
	cd workdir12
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	mafftpath=$type/61mafft_corrected
	selmafftpath=$type/$selection/61mafft_corrected
	alnpath=$type/80concatenated_exon_alignments_corrected
	selalnpath=$type/$selection/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
	echo -e "...with corrected reading frame\n"
else
	mafftpath=$type/60mafft
	selmafftpath=$type/$selection/60mafft
	alnpath=$type/70concatenated_exon_alignments
	selalnpath=$type/$selection/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
	echo -e "\n"
fi

#Check input files
echo -ne "Testing if input data are available..."
if [ -d "$path/$mafftpath" ]; then
	if [ "$(ls -A $path/$mafftpath)" ]; then
		if [ -f "$source/excludelist.txt" ]; then
			echo -e "OK\n"
		else
			echo -e "'excludelist.txt' is missing in 'HybSeqSource'. Exiting...\n"
			rm -d ../workdir12/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/$mafftpath' is empty. Exiting...\n"
		rm -d ../workdir12/ 2>/dev/null
		exit 3
	fi
else
	echo -e "'$path/$mafftpath' is missing. Exiting...\n"
	rm -d ../workdir12/ 2>/dev/null
	exit 3
fi

#Check '$selection'
if [ -z "$selection" ]; then
	echo empty
	echo -e "variable 'selection' is empty. Exiting...\n"
	rm -d ../workdir12/ 2>/dev/null
	exit 3
else
	if ! [[ "$selection" =~ [^a-zA-Z0-9] ]]; then
		echo -e "Modified alignment will be saved to '$type/$selection'"
	else
		echo -e "variable 'selection' contains invalid characters. Exiting...\n"
		rm -d ../workdir12/ 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [ -d "$path/$selmafftpath" ]; then
	echo -e "Directory '$path/$selmafftpath' already exists. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir12/ 2>/dev/null
	exit 3
else
	if [ -d "$path/$selalnpath" ]; then
		echo -e "Directory '$path/$selalnpath' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir12/ 2>/dev/null
		exit 3
	else
		if [[ ! $location == "1" ]]; then
			if [ "$(ls -A ../workdir12)" ]; then
				echo -e "Directory 'workdir12' already exists. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir12/ 2>/dev/null
				exit 3
			fi
		fi
	fi
fi

#Write log
logname=HPM12
echo -e "HybPhyloMaker12: exclude samples" > ${logname}.log
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
for set in data selection cp corrected OUTGROUP; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Check $OUTGROUP against excludelist.txt (if outgroup is to be selected ask for a new outgroup and exit)
#download exclude list (excludelist.txt) from HybSeqSource
cp $source/excludelist.txt .
if [ ! -z $(grep "$OUTGROUP" excludelist.txt) ]; then
	echo -e "Outgroup is in 'excludelist.txt'. Redefine the OUTGROUP or leave it empty. Exiting...\n"
	rm -d ../workdir12/ 2>/dev/null
	exit 3
fi

#Make a new folder for results
mkdir -p $path/$type/$selection
#copy excludelist.txt to $type/$selection
cp excludelist.txt $path/$type/$selection

#Modify 'excludelist.txt'
#change EOLs in excludelist to Unix style
sed -i.bak 's/\x0D$//' excludelist.txt
#remove empty lines from exclude list (very important, otherwise 'grep -v -f' is not working properly!!!)
sed -i.bak2 '/^$/d' excludelist.txt
#add LF at the end of last line in excludelist.txt if missing
sed -i.bak3 '$a\' excludelist.txt
#add '_contigs.fas' or '_contigs_cpDNA.fas' to the end of each line in excludelist.txt (to match names in mafft files)
if [[ $cp =~ "yes" ]]; then
	sed -i.bak4 's/$/_contigs_cpDNA.fas/' excludelist.txt
else
	sed -i.bak4 's/$/_contigs.fas/' excludelist.txt
fi

#Download alignments from '60mafft' folder
if [[ $cp =~ "yes" ]]; then
	cp $path/$mafftpath/*.fasta .
else
	if [[ $corrected =~ "yes" ]]; then
		find $path/$mafftpath -maxdepth 1 -name "*.corrframe" -exec cp -t . {} + #to avoid 'Argument list too long' error
		#remove '.corrframe' (files have now '.mafft' only)
		for filename in $(ls *.corrframe); do
			mv ${filename} ${filename%.*}
		done
	else
		find $path/$mafftpath -maxdepth 1 -name "*.mafft" -exec cp -t . {} + #to avoid 'Argument list too long' error
	fi
	#remove '.mafft' (files have now '.fasta' only)
	for filename in $(ls *.mafft); do
		mv ${filename} ${filename%.*}
	done
fi

#Change leading/tailing '-', i.e. differentiate missing data from gaps
#make a new folder for results
mkdir $path/$selmafftpath
echo -ne "\nRemoving samples from exon alignments..."
ls | grep '.fasta' > listOfMAFFTFiles.txt
for mafftfile in $(cat listOfMAFFTFiles.txt); do
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $mafftfile > tmp && mv tmp $mafftfile
	#Replace leading and tailing '-' by '?'
	perl -pe 's/\G-|-(?=-*$)/?/g' $mafftfile > tmp && mv tmp $mafftfile
	#remove EOL after sample to exclude, i.e. put 'sample to exclude' and its sequence on a single line
	for exclude in $(cat excludelist.txt); do
		perl -pe "s/$exclude\n/$exclude/g" $mafftfile > tmp
		mv tmp $mafftfile
	done
	#get only lines without (-v) samples to exclude
	grep -v -f excludelist.txt $mafftfile > tmp
	if [[ $cp =~ "no" ]]; then
		if [[ $corrected =~ "yes" ]]; then
			#add '.mafft.corrframe' to conform with filenames in 61mafft_corrected folders
			mv tmp ${mafftfile}.mafft.corrframe
			cp ${mafftfile}.mafft.corrframe $path/$selmafftpath
		else
			#add '.mafft' to conform with filenames in 60mafft folders
			mv tmp ${mafftfile}.mafft
			cp ${mafftfile}.mafft $path/$selmafftpath
		fi
	else
		mv tmp ${mafftfile}
		cp ${mafftfile} $path/$selmafftpath
	fi
done


echo -e "finished"

#Concatenate exon alignments
if [[ $cp =~ "no" ]]; then
	echo -ne "\nConcatenating exons to loci..."
	#Copy script
	cp $source/AMAS.py .
	
	if [[ $corrected =~ "yes" ]]; then
		#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft.corrframe to To_align_Assembly_10372_*mafft.corrframe)
		#(all files starting with "To_align_Assembly_10372_" will be merged)
		ls -1 | grep 'corrframe' | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft.corrframe/g > fileNamesForConcat.txt
		#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. CorrectedAssembly_10372
		ls -1 | grep 'corrframe' | cut -d'_' -f4 | sort -u | sed s/^/CorrectedAssembly_/g > fileNamesForSaving.txt
	else
		#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft to To_align_Assembly_10372_*mafft)
		#(all files starting with "To_align_Assembly_10372_" will be merged)
		ls -1 | grep 'mafft' | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft/g > fileNamesForConcat.txt
		#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. Assembly_10372
		ls -1 | grep 'mafft' | cut -d'_' -f4 | sort -u | sed s/^/Assembly_/g > fileNamesForSaving.txt
	fi
	#Combine both files (make single file with two columns)
	paste fileNamesForConcat.txt fileNamesForSaving.txt > fileForLoop.txt
	
	if [[ $location == "1" ]]; then
		#Add necessary module
		module add python-3.4.1-gcc
		#Make a new folder for results
		mkdir $path/$selalnpath
	elif [[ $location == "2" ]]; then
		#Add necessary module
		module unload bioinformatics/anaconda3 #unload possible previously loaded python3
		module load bioinformatics/anaconda3/2.3.0 #python3
		#Make a new folder for results
		mkdir $path/$selalnpath
	else
		mkdir $path/$selalnpath
	fi
	#Concatenate the exon alignments (values from first and second column of fileForLoop.txt are assigned to variable 'a' and 'b', respectively),
	#transform fasta to phylip format, copy results from scratch to home
	cat fileForLoop.txt | while read -r a b
	do
		#concatenate exons from the same locus
		python3 AMAS.py concat -i $a -f fasta -d dna -u fasta -t ${b}.fasta >/dev/null
		python3 AMAS.py concat -i $a -f fasta -d dna -u phylip -t ${b}.phylip >/dev/null
		#modify and rename partitions.txt
		sed -i.bak -e 's/^/DNA, /g' -e 's/_To_align//g' partitions.txt
		mv partitions.txt ${b}.part
		if [[ $location == "1" ]]; then
			cp $b.* $path/$selalnpath
		else
			cp $b.* $path/$selalnpath
		fi
	done
	echo -e "finished"
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/$selection

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir12
fi

echo -e "\nHybPhyloMaker12 finished...\n"

