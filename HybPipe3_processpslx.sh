#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=8gb
#PBS -N HybPipe3_process_pslx
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe3_process_pslx
#$ -o HybPipe3_process_pslx.log


# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                        Script 03 - Process pslx files                        *
# *                                   v.1.0.3                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#Input: pslx files named genus-species_code_contigs.fas.pslx in the folder $otherpslx

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo "Metacentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	#Add necessary modules
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
	mkdir workdir03
	cd workdir03
	#Add necessary modules
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
	mkdir workdir03
	cd workdir03
fi

#-----------------------COMBINATION OF SEQUENCES OF THE EXONS OF EACH ACCESSION-----------------------
#Copy script and reference
cp -r $source/assembled_exons_to_fastas.py .
cp -r $source/$probes .
chmod +x assembled_exons_to_fastas.py
if [[ $location == "1" ]]; then
	#Add necessary modules
	module add python-2.7.6-gcc
	module add python-2.7.6-intel
fi
#Copy other pslx files to combine
cp -r $otherpslxpath/* .
#Make a list of all pslx files
ls *.pslx > listOfPslxFiles.txt
#Run script that creates folder 'contigsMatchLoci' which contains files for all exons
python assembled_exons_to_fastas.py -l listOfPslxFiles.txt -f $probes -d contigsMatchLoci
echo "Finished combination of sequences..."

#-----------------------ALIGNING FASTA FILES (ALL EXONS FOR ALL SPECIES) USING MAFFT-----------------------
#Enter directory with fasta files
cd contigsMatchLoci
if [[ $location == "1" ]]; then
	#Add necessary module
	module add mafft-7.029
	module add parallel
elif [[ $location == "2" ]]; then
	module load bioinformatics/mafft/7.221
	module load tools/gnuparallel/20160422
fi
#Make a list of all fasta files
ls *.fasta > listOfFastaFiles.txt
#Make a new folder for results
if [[ $location == "1" ]]; then
	mkdir $path/60mafft
else
	mkdir ../$path/60mafft
fi
#A loop/parallelization to process all samples in folders named as specified in listOfFastaFiles.txt
if [ "$parallelmafft" = "yes" ]; then
	cat listOfFastaFiles.txt | parallel -j 8 --max-procs 8 'mafft --auto {} > {}.mafft'
	if [[ $location == "1" ]]; then
		cp *.mafft $path/60mafft
	else
		cp *.mafft ../$path/60mafft
	fi
else
	for fastafile in $(cat listOfFastaFiles.txt)
	do
		mafft --auto $fastafile > $fastafile.mafft
		if [[ $location == "1" ]]; then
			cp $fastafile.mafft $path/60mafft
		else
			cp $fastafile.mafft ../$path/60mafft
		fi
	done
fi
echo "Finished MAFFT alignment..."

#-----------------------CHANGE LEADING AND TAILING '-' TO '?'-----------------------
#i.e. differentiate missing data from gaps
ls *.mafft > listOfMAFFTFiles.txt
for mafftfile in $(cat listOfMAFFTFiles.txt)
do
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $mafftfile > tmp && mv tmp $mafftfile
	#Replace leading and tailing '-' by '?'
	sed -i -e ':a;s/^\(-*\)-/\1?/;ta' -e ':b;s/-\(-*\)$/?\1/;tb' $mafftfile
done
echo "Finished changing leading/tailing '-' ..."

#-----------------------CONCATENATE THE EXON ALIGNMENTS-----------------------
#Copy script
if [[ $location == "1" ]]; then
	cp -r $source/catfasta2phyml.pl .
else
	cp -r ../$source/catfasta2phyml.pl .
fi
#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft to To_align_Assembly_10372_*mafft)
#(all files starting with "To_align_Assembly_10372_" will be merged)
ls -1 *.mafft | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft/g > fileNamesForConcat.txt
#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. Assembly_10372
ls -1 *.mafft | cut -d'_' -f4 | sort -u | sed s/^/Assembly_/g > fileNamesForSaving.txt
#Combine both files (make single file with two columns)
paste fileNamesForConcat.txt fileNamesForSaving.txt > fileForLoop.txt

if [[ $location == "1" ]]; then
	#Add necessary module
	module add perl-5.10.1
	#Make a new folder for results
	mkdir $path/70concatenated_exon_alignments
else
	mkdir ../$path/70concatenated_exon_alignments
fi
#Concatenate the exon alignments (values from first and second column of fileForLoop.txt are assigned to variable 'a' and 'b', respectively),
#transform fasta to phylip format, copy results from scratch to home
cat fileForLoop.txt | while read -r a b
do
	perl catfasta2phyml.pl -f $a > $b.fasta
	perl catfasta2phyml.pl $b.fasta > $b.phylip
	if [[ $location == "1" ]]; then
		cp $b.* $path/70concatenated_exon_alignments
	else
		cp $b.* ../$path/70concatenated_exon_alignments
	fi
done
echo "Finished concatenation of exons..."
#Move back to workdir
cd ..

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir03
fi
