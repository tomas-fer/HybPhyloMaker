#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:0:0
#PBS -l select=1:ncpus=2:mem=1gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker4a_selectNonHet
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker4a2_selectNonHet
#$ -o HybPhyloMaker4a2_selectNonHet.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                 Script 04a2 - Select low-heterozygosity exons                *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

#Input:
#- exons/21_mapped_{mappingmethod}
#- exons/60mafft

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker4a2 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	otherpslxcppath=/storage/$server/home/$LOGNAME/$otherpslxcp
	#Add necessary modules
	#module add python-2.7.6-gcc
	#module add python-2.7.6-intel
	module add mafft-7.029
	module add parallel
	module add perl-5.10.1
	module add bwa-0.7.15
	module add samtools-1.3
	module add bcftools-1.3.1
	module add debian8-compat
	#available processors
	procavail=`expr $TORQUE_RESC_TOTAL_PROCS - 2`
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker4a2 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	otherpslxcppath=../$otherpslxcp
	#Make and enter work directory
	mkdir -p workdir04a2
	cd workdir04a2
	#Add necessary modules
	module load bioinformatics/mafft/7.221
	module load tools/gnuparallel/20160422
else
	echo -e "\nHybPhyloMaker4a2 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	otherpslxcppath=../$otherpslxcp
	#Make and enter work directory
	mkdir -p workdir04a2
	cd workdir04a2
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
else
	echo -e "Working with exons\n"
	type="exons"
fi

#Check necessary files
echo -ne "Testing if input data are available..."
if [[ $cp =~ "yes" ]]; then
	if [ -f "$source/$cpDNACDS" ]; then
		if [ -d "$otherpslxcppath" ]; then
			if [ "$(ls -A $otherpslxcppath)" ]; then
				echo -e "OK\n"
			else
				echo -e "'$otherpslxcppath' is empty. Move desired *.pslx files into it.\nExiting...\n"
				rm -d ../workdir04a2/ 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$otherpslxcppath' does not exists. Create this directory and move desired *.pslx files into it.\nExiting...\n"
			rm -d ../workdir04a2/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$cpDNACDS' is missing in 'HybSeqSource'. Exiting...\n"
		rm -d ../workdir04a2/ 2>/dev/null
		exit 3
	fi
else
	if [ -f "$source/$probes" ]; then
		if [ -d "$otherpslxpath" ]; then
			if [ "$(ls -A $otherpslxpath)" ]; then
				echo -e "OK\n"
			else
				echo -e "'$otherpslxpath' is empty. Move desired *.pslx files into it.\nExiting...\n"
				rm -d ../workdir04a2/ 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$otherpslxpath' does not exists. Create this directory and move desired *.pslx files into it.\nExiting...\n"
			rm -d ../workdir04a2/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$probes' is missing in 'HybSeqSource'. Exiting...\n"
		rm -d ../workdir04a2/ 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [ ! -d "$path/$type/60mafft" ]; then
	echo -e "Directory '$path/$type/60mafft' does not exists. Run first HybPhyloMaker4a. Exiting...\n"
	rm -d ../workdir04a2/ 2>/dev/null
	exit 3
else
	if [ -d "$path/$type/70concatenated_exon_alignments_NoHet" ]; then
		if [[ $nohetcalculation =~ "yes" ]]; then
			echo -e "Directory '$path/$type/70concatenated_exon_alignments_NoHet' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir04a2/ 2>/dev/null
			exit 3
		elif [ -d "$path/$type/70concatenated_exon_alignments_NoHet/${maxhet}" ]; then
			echo -e "Directory '$path/$type/70concatenated_exon_alignments_NoHet/${maxhet}' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir04a2/ 2>/dev/null
			exit 3
		fi
	else
		if [[ ! $location == "1" ]]; then
			if [ "$(ls -A ../workdir04a2)" ]; then
				echo -e "Directory 'workdir04a2' already exists. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir04a2/ 2>/dev/null
				exit 3
			fi
		fi
	fi
fi

#Test if summary of heterozygous sites already exists
if [[ $nohetcalculation =~ "no" ]] && [ ! -f "$path/$type/70concatenated_exon_alignments_NoHet/summary/summary.txt" ]; then
	echo -e "Summary file does not exists. Run this script with 'nohetcalculation=yes' in order to create it. Exiting...\n"
	rm -d ../workdir04a2/ 2>/dev/null
	exit 3
fi

#Write log
logname=HPM4a2
echo -e "HybPhyloMaker4a2: select low-heterozygosity exons" > ${logname}.log
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
for set in data cp probes cpDNACDS otherpslx otherpslxcp maxhet; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	cp $source/$cpDNACDS .
else
	cp $source/$probes .
fi

# Make new folder for results
mkdir -p $path/$type/70concatenated_exon_alignments_NoHet
mkdir -p $path/$type/70concatenated_exon_alignments_NoHet/summary
mkdir -p $path/$type/70concatenated_exon_alignments_NoHet/${maxhet}

if [[ $nohetcalculation =~ "yes" ]]; then
	#--------------VARIANT CALLING AND CALCULATION OF NUMBER OF HETEROZYGOUS SITES PER EXON---------------
	# Cut filename before first '.', i.e. remove suffix - does not work if there are other dots in reference file name
	if [[ $cp =~ "yes" ]]; then
		name=`ls $cpDNACDS | cut -d'.' -f 1`
	else
		name=`ls $probes | cut -d'.' -f 1`
	fi
	reference=${name}_with${nrns}Ns_beginend.fas
	cp $source/$reference .
	header=$(grep ">" ${name}_with${nrns}Ns_beginend.fas | sed 's/>//')
	bedfile=${header}_exon_positions.bed
	cp $source/$bedfile .
	# Get exon positions
	awk '{ print $2 "\t" $3}' $bedfile > positions.txt
	# Get exon names
	grep ">" $probes > summary.txt
	paste summary.txt positions.txt > tmp && mv tmp summary.txt
	#add as first line
	sed -i "1i exon\tFrom\tTo" summary.txt
	sed -i 's/>//' summary.txt
	cp summary.txt commoncolumns.txt
	# Copy BAM files
	cp $path/$type/21mapped_${mappingmethod}/*.bam .
	#indexing reference
	bwa index $reference 2>/dev/null
	# Sorting, indexing, variant calling, extraxcting heteroyzgous positions, calculating nr hetero per exon
	for file in $(ls *.bam | cut -d'.' -f1); do
		echo -e "$file"
		echo -e "...sorting"
		samtools sort ${file}.bam -o ${file}_sorted.bam
		rm ${file}.bam
		echo -e "...indexing"
		samtools index ${file}_sorted.bam
		#samtools view -h -o ${file}.sam ${file}_sorted.bam
		echo -e "...pileup"
		samtools mpileup -E -uf $reference ${file}_sorted.bam > ${file}.pileup 2>/dev/null
		echo -e "...vcf"
		bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v -o ${file}.vcf ${file}.pileup 2>/dev/null
		#bcftools view -v -m 0.5 ${file}.pileup > ${file}.vcf
		rm ${file}.pileup
		#grep heterozygous
		echo -e "...grep hetero"
		grep -e 0/1 -e 1/2 ${file}.vcf > ${file}_hetero.txt
		#only SNPs with mapping quality higher than 36 (mapping quality in 6th column)
		awk '{if ($6 >= 36) print $0}' ${file}_hetero.txt > ${file}_heteroqual.txt
		echo -e "...calculate nr het per exon\n"
		echo -e "$file" > ${file}_nrhet.txt
		cat $bedfile | while read a b c; do
			nrhet=$(awk -v low=$b -v high=$c '{ if (($2 >= low) && ($2 <= high)) print $0 }' ${file}_heteroqual.txt | wc -l)
			echo -e "$nrhet" >> ${file}_nrhet.txt
		done
		paste commoncolumns.txt ${file}_nrhet.txt > ${file}_nrhet_per_exon.txt
		cat ${file}_nrhet.txt | cut -f3 > ${file}_nrhetonly.txt
		paste summary.txt ${file}_nrhetonly.txt > tmp && mv tmp summary.txt
		rm ${file}_nrhetonly.txt ${file}_nrhet.txt
		sed -i 's/>//' ${file}_nrhet_per_exon.txt
		cp ${file}_nrhet_per_exon.txt $path/$type/70concatenated_exon_alignments_NoHet/summary
		#only exons with maximum specified heterozygosity
		awk -v val=${maxhet} '$4 <= val' ${file}_nrhet_per_exon.txt | cut -d'_' -f1,2 | sort | uniq > ${file}_ExonsWithoutHet${maxhet}.txt
		cp ${file}_ExonsWithoutHet${maxhet}.txt $path/$type/70concatenated_exon_alignments_NoHet/summary
	done
	
	#extract first line (header) and add two columns to the end 'total' and 'max'
	head -n1 summary.txt | awk '{ print $0 "\t" "total" "\t" "max" }' > headerSummary.txt
	sed -i '1d' summary.txt # delete first line
	# Compute sum from numbers on each line (omitting first three columns with assembly name, start and stop positions), i.e. sum of heterozygous sites over all samples,
	# and add this as a last value on each line (last sed command replaces double spaces by TABs)
	awk '
	BEGIN {FS=OFS=" "}
	{
	sum=0; n=0
	for(i=4;i<=NF;i++)
		 {sum+=$i; ++n}
		 print $0," "sum
		 sum=0; n=0
	}' summary.txt | sed "s/  /\t/" > tmp && mv tmp summary.txt
	#Add highest number of heterozygous sites per locus over samples as the last column in 'summary.txt'
	awk '{max=$4;for(i=4;i<=NF-1;i++){if($i > max) max = $i}print $0 "\t"max}' summary.txt > tmp && mv tmp summary.txt
	cat headerSummary.txt summary.txt > tmp && mv tmp summary.txt
	cp summary.txt $path/$type/70concatenated_exon_alignments_NoHet/summary
else
	cp $path/$type/70concatenated_exon_alignments_NoHet/summary/summary.txt .
fi

#-----------CREATE A LIST OF EXONS WITH LIMITED NUMBER OF NON-HETEROZYGOUS SITES OVER ALL SAMPLES-----------
awk -v val=${maxhet} '$NF <= val { print $1 }' summary.txt | sort | uniq > exonsWithoutHet${maxhet}.txt
cp exonsWithoutHet${maxhet}.txt $path/$type/70concatenated_exon_alignments_NoHet/${maxhet}

#------------------------------COPY MAFFT FILES-------------------------------
for i in $(cat exonsWithoutHet${maxhet}.txt); do
	cp $path/$type/60mafft/To_align_${i}.fasta.mafft .
done

#-----------------------CHANGE LEADING AND TAILING '-' TO '?'-----------------------
#i.e. differentiate missing data from gaps
echo -ne "Changing leading/tailing '-' in alignments..."
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
	#sed -i.bak -e ':a;s/^\(-*\)-/\1?/;ta' -e ':b;s/-\(-*\)$/?\1/;tb' $mafftfile
	perl -pe 's/\G-|-(?=-*$)/?/g' $mafftfile > tmp && mv tmp $mafftfile
	if [[ $cp =~ "yes" ]]; then
		cp $mafftfile $path/$type/60mafft
	fi
done
echo -e "finished"

#-----------------------CONCATENATE THE EXON ALIGNMENTS-----------------------
if [[ $cp =~ "no" ]]; then
	echo -ne "\nConcatenating exons to loci..."
	#Copy script
	#cp -r $source/catfasta2phyml.pl .
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
		#module add perl-5.10.1
	elif [[ $location == "2" ]]; then
		#Add necessary module
		module unload bioinformatics/anaconda3 #unload possible previously loaded python3
		module load bioinformatics/anaconda3/2.3.0 #python3
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
		#perl catfasta2phyml.pl -f $a > $b.fasta
		#perl catfasta2phyml.pl $b.fasta > $b.phylip
		cp $b.* $path/$type/70concatenated_exon_alignments_NoHet/${maxhet}
	done
	echo -e "finished"
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/70concatenated_exon_alignments_NoHet

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir04a2
fi

echo -e "\nScript HybPhyloMaker4a2 finished...\n"
