#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker4a3_exonTrimming
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 2
#$ -q mThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker4a3_removeEmptyPositions
#$ -o HybPhyloMaker4a3_removeEmptyPositions.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                          Script 4a3 - TrimmingExons                          *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker4a3 is running on MetaCentrum...\n"
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
	module add perl-5.10.1
	module add trimal-1.4
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker4a3 is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir01
	cd workdir01
	#Add necessary modules
	module load bioinformatics/trimal/1.4
else
	echo -e "\nHybPhyloMaker4a3 is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir04a3
	cd workdir04a3
fi

#Write log
logname=HPM4a3
echo -e "HybPhyloMaker4a3: exon trimming" > ${logname}.log
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
for set in data cp gappyout noallgaps; do
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

if [[ $noallgaps =~ "no" && $gappyout =~ "no" ]] ; then
	echo -e "Nothing to do. Exiting..."
	echo -e "Set 'noallgaps' and/or 'gappyout' to 'yes' and rerun.\n"
	rm -d ../workdir04a3/ 2>/dev/null
	exit 3
fi

#Copy assemblies to scratch/workdir
if [[ $cp =~ "yes" ]]; then
	cp $path/$type/60mafft/*.fasta .
	ls *.fasta > listOfMAFFTFiles.txt
else
	cp $path/$type/60mafft/*.mafft .
	ls | grep '.mafft' > listOfMAFFTFiles.txt
fi

#make folder for results
mkdir -p $path/$type/60mafft/trimmed
mkdir -p $path/$type/60mafft/trimmed/html

#PROCESS MAFFT FILES
echo -e "\nTrimming alignments..."
numbermafft=$(cat listOfMAFFTFiles.txt | wc -l)
calculating=0
for mafft in $(cat listOfMAFFTFiles.txt); do
	calculating=$((calculating + 1))
	echo -e "\t${mafft} ($calculating out of $numbermafft)"
	echo ${mafft} >> warnings.txt
	#trimming (keep sequences even if only composed by gaps)
	if [[ $noallgaps =~ "yes" ]]; then
		trimal -in ${mafft} -out ${mafft}.1 -htmlout ${mafft}.noallgaps.html -noallgaps -keepseqs >/dev/null 2>> warnings.txt #remove positions with gaps only
		mv ${mafft}.1 ${mafft}
	fi
	if [[ $gappyout =~ "yes" ]]; then
		trimal -in ${mafft} -out ${mafft}.2 -htmlout ${mafft}.gappyout.html -gappyout -keepseqs >/dev/null 2>> warnings.txt #remove gappy positions
		mv ${mafft}.2 ${mafft}
	fi
	#change '-' to 'Q'
	sed -i '/^>/!s/-/Q/g' ${mafft}
	#change 'n' to '-' (to allow trimming on Ns)
	sed -i '/^>/!s/n/-/g' ${mafft}
	#second trimming (keep sequences even if only composed by gaps)
	if [[ $noallgaps =~ "yes" ]]; then
		trimal -in ${mafft} -out ${mafft}.3 -htmlout ${mafft}.noallgapsN.html -noallgaps -keepseqs >/dev/null 2>> warnings.txt #remove positions with gaps only
		mv ${mafft}.3 ${mafft}
	fi
	if [[ $gappyout =~ "yes" ]]; then
		trimal -in ${mafft} -out ${mafft}.4 -htmlout ${mafft}.gappyoutN.html -gappyout -keepseqs >/dev/null 2>> warnings.txt #remove gappy positions
		mv ${mafft}.4 ${mafft}
	fi
	#change '-' back to 'n'
	sed -i '/^>/!s/-/n/g' ${mafft}
	#change 'Q' back to '-'
	sed -i '/^>/!s/Q/-/g' ${mafft}
	#removes everything after a gap (i.e., also sequence length introduced by trimAl)
	sed -i 's/ .*//' ${mafft}
	cp ${mafft} $path/$type/60mafft/trimmed
	cp ${mafft}.*.html $path/$type/60mafft/trimmed/html
	rm ${mafft}.*.html
done
echo -e "finished"

#Remove empty files locally
find . -name '*.mafft' -type f -empty -delete
#Remove empty file in home dir
find $path/$type/60mafft/trimmed/ -name '*.mafft' -type f -empty -delete

#-----------------------CHANGE LEADING AND TAILING '-' TO '?'-----------------------
#i.e. differentiate missing data from gaps
echo -ne "\nChanging leading/tailing '-' in alignments..."
for mafftfile in $(cat listOfMAFFTFiles.txt); do
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $mafftfile > tmp && mv tmp $mafftfile
	#Replace leading and tailing '-' by '?'
	#sed -i.bak -e ':a;s/^\(-*\)-/\1?/;ta' -e ':b;s/-\(-*\)$/?\1/;tb' $mafftfile
	perl -pe 's/\G-|-(?=-*$)/?/g' $mafftfile > tmp && mv tmp $mafftfile
	if [[ $cp =~ "yes" ]]; then
		cp $mafftfile $path/$type/60mafft/trimmed
	fi
done
echo -e "finished"

#-----------------------CONCATENATE THE EXON ALIGNMENTS-----------------------
if [[ $cp =~ "no" ]]; then
	echo -ne "\nConcatenating exons to loci..."
	#Copy script
	cp -r $source/AMAS.py .
	#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft to To_align_Assembly_10372_*mafft)
	#(all files starting with "To_align_Assembly_10372_" will be merged)
	ls -1 | grep 'mafft' | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft/g > fileNamesForConcat.txt
	#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. Assembly_10372
	ls -1 | grep 'mafft' | cut -d'_' -f4 | sort -u | sed s/^/Assembly_/g > fileNamesForSaving.txt
	#Combine both files (make single file with two columns)
	paste fileNamesForConcat.txt fileNamesForSaving.txt > fileForLoop.txt

	if [[ $location == "1" ]]; then
		#Add necessary module
		module add python-3.4.1-intel
		#Make a new folder for results
		mkdir $path/$type/70concatenated_exon_alignments/trimmed
	elif [[ $location == "2" ]]; then
		#Add necessary module
		module unload bioinformatics/anaconda3 #unload possible previously loaded python3
		module load bioinformatics/anaconda3/2.3.0 #python3
		#Make a new folder for results
		mkdir $path/$type/70concatenated_exon_alignments/trimmed
	else
		mkdir $path/$type/70concatenated_exon_alignments/trimmed
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
		cp $b.* $path/$type/70concatenated_exon_alignments/trimmed
	done
	echo -e "finished"
fi

#Copy warnings to mafft
cp warnings.txt $path/$type/60mafft/trimmed

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/70concatenated_exon_alignments/trimmed

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir04a3
fi

echo -e "\nScript HybPhyloMaker4a3 finished...\n"

