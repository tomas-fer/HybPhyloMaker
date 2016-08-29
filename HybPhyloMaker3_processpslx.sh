#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=8gb
#PBS -N HybPhyloMaker3_process_pslx
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker3_process_pslx
#$ -o HybPhyloMaker3_process_pslx.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                        Script 03 - Process pslx files                        *
# *                                   v.1.1.2                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#Input: pslx files named genus-species_code_contigs.fas.pslx in the folder $otherpslx

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker3 is running on MetaCentrum...\n"
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
	otherpslxcppath=/storage/$server/home/$LOGNAME/$otherpslxcp
	#Add necessary modules
	module add python-2.7.6-gcc
	module add python-2.7.6-intel
	module add mafft-7.029
	module add parallel
	module add perl-5.10.1
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker3 is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	otherpslxcppath=../$otherpslxcp
	#Make and enter work directory
	mkdir workdir03
	cd workdir03
	#Add necessary modules
	module load bioinformatics/mafft/7.221
	module load tools/gnuparallel/20160422
else
	echo -e "\nHybPhyloMaker3 is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	otherpslxcppath=../$otherpslxcp
	#Make and enter work directory
	mkdir workdir03
	cd workdir03
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="_cp"
else
	echo -e "Working with exons\n"
	type=""
fi

#-----------------------COMBINATION OF SEQUENCES OF THE EXONS OF EACH ACCESSION-----------------------
echo -ne "Combining sequences..."
#Copy script and reference
cp -r $source/assembled_exons_to_fastas.py .
if [[ $cp =~ "yes" ]]; then
	cp -r $source/$cpDNACDS .
else
	cp -r $source/$probes .
fi
chmod +x assembled_exons_to_fastas.py
# if [[ $location == "1" ]]; then
	# #Add necessary modules
	# module add python-2.7.6-gcc
	# module add python-2.7.6-intel
# fi

if [[ $cp =~ "yes" ]]; then
	#Copy other cpDNA pslx files to combine
	cp -r $otherpslxcppath/* .
else
	#Copy other pslx files to combine
	cp -r $otherpslxpath/* .
fi

#Make a list of all pslx files
ls *.pslx > listOfPslxFiles.txt
#Run script that creates folder 'contigsMatchLoci' which contains files for all exons
if [[ $cp =~ "yes" ]]; then
	python assembled_exons_to_fastas.py -l listOfPslxFiles.txt -f $cpDNACDS -d contigsMatchLoci
else
	python assembled_exons_to_fastas.py -l listOfPslxFiles.txt -f $probes -d contigsMatchLoci
fi
echo -e "finished"

#-----------------------ALIGNING FASTA FILES (ALL EXONS FOR ALL SPECIES) USING MAFFT-----------------------
#Enter directory with fasta files
cd contigsMatchLoci
# if [[ $location == "1" ]]; then
	# #Add necessary module
	# module add mafft-7.029
	# module add parallel
# elif [[ $location == "2" ]]; then
	# module load bioinformatics/mafft/7.221
	# module load tools/gnuparallel/20160422
# fi
#Make a list of all fasta files
ls *.fasta > listOfFastaFiles.txt
#Make a new folder for results
if [[ $location == "1" ]]; then
	mkdir $path/60mafft${type}
else
	mkdir ../$path/60mafft${type}
fi
#A loop/parallelization to process all samples in folders named as specified in listOfFastaFiles.txt
if [[ $cp =~ "yes" ]]; then
	if [ "$parallelmafft" = "yes" ]; then
		echo -e "\nAligning exons using MAFFT with GNU parallel on $numbcores cores..."
		if [[ $location == "1" ]]; then
			cat listOfFastaFiles.txt | parallel -j $TORQUE_RESC_TOTAL_PROCS 'mafft --auto {} > {}.mafft'
		elif [[ $location == "2" ]]; then
			cat listOfFastaFiles.txt | parallel -j $NSLOTS --max-procs $NSLOTS 'mafft --auto {} > {}.mafft'
		else
			echo -e "Check progress of parallel jobs below\n"
			cat listOfFastaFiles.txt | parallel --eta -j $numbcores 'mafft --auto {} > {}.mafft 2>/dev/null'
		fi
		# if [ ! $LOGNAME == "" ]; then
			# cp *.mafft $path/60mafft_cp
		# else
			# cp *.mafft ../$path/60mafft_cp
		# fi
	else
		echo -e "\nAligning exons using MAFFT one by one..."
		for fastafile in $(cat listOfFastaFiles.txt); do
			echo $fastafile
			mafft --auto $fastafile > $fastafile.mafft 2>/dev/null
			# if [ ! $LOGNAME == "" ]; then
				# cp $fastafile.mafft $path/60mafft_cp
			# else
				# cp $fastafile.mafft ../$path/60mafft_cp
			# fi
		done
	fi
	echo -e "\nFinished MAFFT alignment..."
	#Remove all fasta files (to be able to work only with aligned fasta files (renamed mafft files - see next step)
	rm *.fasta
	#Rename mafft files (to Assembly_NAME.fasta) and if the NAME already exists add 'x2' (to get Assembly_NAMEx2.fasta)
	for i in *.mafft; do
		if [ -f Assembly_`basename "$i" .mafft | cut -f5 -d "_"` ]; then
			#if the filename after renaming is found, add 'x2' after the gene name
			mv "$i" Assembly_`basename "$i" .mafft | cut -f5 -d "_" | sed 's/\./x2\./'`
		else
			#change the name
			mv "$i" Assembly_`basename "$i" .mafft | cut -f5 -d "_"`
		fi
	done
else
	if [ "$parallelmafft" = "yes" ]; then
		echo -e "\nAligning exons using MAFFT with GNU parallel on $numbcores cores..."
		if [[ $location == "1" ]]; then
			cat listOfFastaFiles.txt | parallel -j $TORQUE_RESC_TOTAL_PROCS 'mafft --auto {} > {}.mafft'
		elif [[ $location == "2" ]]; then
			cat listOfFastaFiles.txt | parallel -j $NSLOTS --max-procs $NSLOTS 'mafft --auto {} > {}.mafft'
		else
			echo -e "Check progress of parallel jobs below\n"
			cat listOfFastaFiles.txt | parallel --eta -j $numbcores 'mafft --auto {} > {}.mafft 2>/dev/null'
		fi
		if [[ $location == "1" ]]; then
			cp *.mafft $path/60mafft
		else
			find . -name "*.mafft" -exec cp -t ../$path/60mafft/ {} +
		fi
	else
		echo -e "\nAligning exons using MAFFT one by one..."
		for fastafile in $(cat listOfFastaFiles.txt); do
			echo $fastafile
			mafft --auto $fastafile > $fastafile.mafft 2>/dev/null
			if [[ $location == "1" ]]; then
				cp $fastafile.mafft $path/60mafft
			else
				cp $fastafile.mafft ../$path/60mafft
			fi
		done
	fi
	echo -e "\nFinished MAFFT alignment..."
fi
#-----------------------CHANGE LEADING AND TAILING '-' TO '?'-----------------------
#i.e. differentiate missing data from gaps
echo -ne "\nChanging leading/tailing '-' in alignments..."
if [[ $cp =~ "yes" ]]; then
	ls *.fasta > listOfMAFFTFiles.txt
else
	ls | grep '.mafft' > listOfMAFFTFiles.txt
fi
for mafftfile in $(cat listOfMAFFTFiles.txt)
do
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $mafftfile > tmp && mv tmp $mafftfile
	#Replace leading and tailing '-' by '?'
	sed -i.bak -e ':a;s/^\(-*\)-/\1?/;ta' -e ':b;s/-\(-*\)$/?\1/;tb' $mafftfile
	if [[ $cp =~ "yes" ]]; then
		cp $mafftfile ../$path/60mafft_cp
	fi
done
echo -e "finished"

#-----------------------CONCATENATE THE EXON ALIGNMENTS-----------------------
if [[ $cp =~ "no" ]]; then
	echo -ne "\nConcatenating exons to loci..."
	#Copy script
	if [[ $location == "1" ]]; then
		cp -r $source/catfasta2phyml.pl .
	else
		cp -r ../$source/catfasta2phyml.pl .
	fi
	#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft to To_align_Assembly_10372_*mafft)
	#(all files starting with "To_align_Assembly_10372_" will be merged)
	ls -1 | grep 'mafft' | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft/g > fileNamesForConcat.txt
	#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. Assembly_10372
	ls -1 | grep 'mafft' | cut -d'_' -f4 | sort -u | sed s/^/Assembly_/g > fileNamesForSaving.txt
	#Combine both files (make single file with two columns)
	paste fileNamesForConcat.txt fileNamesForSaving.txt > fileForLoop.txt

	if [[ $location == "1" ]]; then
		#Add necessary module
		#module add perl-5.10.1
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
	echo -e "finished"
fi

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

echo -e "\nScript HybPhyloMaker3 finished...\n"
