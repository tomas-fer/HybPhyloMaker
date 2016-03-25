#!/bin/bash
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=4gb
#PBS -l scratch=1gb
#PBS -N HybPipe4_missingdata_handling
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                      Script 04 - Missing data handling                       *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Remove species with given percentage (or more) of missing data (i.e., 'N', 'n' or '?', but NOT '-') from multiple FASTA file
#Removal is done for all *.fasta files in the folder
#All fasta files must include all samples in the same order!!!
#Headers must not include spaces!!!
#Creates following files:
#(1) modified FASTA files without samples with percentage of missing data higher than specified value
#(2) table with headers and percentage of missing data
#(3) list of assemblies with average percentage of missing data bellow specified values and with no completely missing data per assembly
#(4) table of summary characteristics (length, %missing, number of variable and parsimony informative sites, GC content...) for each selected locus
#CUT = delete all samples with this or higher percentage of N

#Complete path and set configuration for selected location
if [ ! $LOGNAME == "" ]; then
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
	#Add necessary modules
	module add python-3.4.1-gcc
else
	echo "Local..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir04
	cd workdir04
fi

#Copy data folder to scratch
cp $path/70concatenated_exon_alignments/*.fasta .

#Make a list of all fasta files
ls *.fasta | cut -d"." -f1 > fileForDeletePercentage.txt
#Make new dir for results
mkdir $path/71selected${CUT}
mkdir $path/71selected${CUT}/deleted_above${CUT}
for file in $(cat fileForDeletePercentage.txt)
do
	#Delete empty lines (to be on the safe side...)
	sed -i '/^$/d' $file.fasta
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file.fasta > tmp && mv tmp $file.fasta
	#Replace newline with ' ' if line starts with '>' (i.e., merge headers with data into single line separated by space)
	sed -i '/^>/{N; s/\n/ /;}' $file.fasta
	
	#Create file for storage headers and percentages of N
	touch ${file}_${CUT}percN.fas
	echo -e "species\t$file" >> ${file}_${CUT}percN.fas
	#Create file for accessions with N below threshold (retained accessions)
	touch ${file}_modif${CUT}.fas
	#Loop for reading file line by line
	cat $file.fasta | while read line
	do
		#Retain only part after the first space, i.e., only sequence
		sequence=$(sed 's/[^ ]* //' <<< $line)
		#Cut first part until space, i.e. header
		header=$(cut -f1 -d" " <<< $line | sed 's/>//' | sed 's/.fas//')
		#Put only 'n', 'N' and '?' from $sequence to $x
		x="${sequence//[^nN?]}"
		#Count number of 'n', 'N' and '?' (i.e., length of $x)
		count="${#x}"
		#Put only '-' (i.e., indels) from $sequence to $x2
		x2="${sequence//[^-]}"
		#Count number of '-' (i.e., length of $x2)
		indel="${#x2}"
		#Count length of $sequence
		a=$sequence
		y=${#a}
		#If $count (missing data) plus $indel (number of '-') equals $y (total sequence length) then $count = $y (i.e., there data completely missing)
		totalmissing=`expr $count + $indel`
		if [ "$totalmissing" -eq "$y" ]; then
			count=$y
		fi
		#Calculate percentage of missing data with 2 decimal points (for printing to the table)
		percentageprint=$(echo -e "scale=2;($count * 100) / $y" | bc)
		#Calculate percentage of missing data as integer (for evaluation in bash)
		percentage=$(echo -e "scale=0;($count * 100) / $y" | bc)
		echo -e "$header\t$percentageprint" >> $file\_$CUT\percN.fas
		if [ "$percentage" -lt "$CUT" ]; then
			echo "$line" >> ${file}_modif${CUT}.fas
		else
			echo "Deleting $header with $percentageprint of missing data"
		fi
	done
	#Replace ' ' by newline character
	sed -i '/^>/{s/ /\n/}' $file.fasta
	sed -i '/^>/{s/ /\n/}' ${file}_modif${CUT}.fas
	#Copy results home
	cp ${file}_${CUT}percN.fas $path/71selected${CUT}/deleted_above${CUT}
	cp ${file}_modif${CUT}.fas $path/71selected${CUT}/deleted_above${CUT}
	echo "$file processed"
done

#Prepare header file
file=$(cat fileForDeletePercentage.txt | head -n 1)
cat $file.fasta | sed '/^>/{N; s/\n/ /;}' | cut -f1 -d" " | sed 's/>//' | sed 's/_contigs//' | sed 's/.fas//' > headers.txt

#Make a table with percentages of N in each accession and file
awk '{_[FNR]=(_[FNR] OFS $2)}END{for (i=1; i<=FNR; i++) {sub(/^ /,"",_[i]); print _[i]}}' *percN.fas > missing_percentage_overview.txt
#Add word 'species' as a 1st line in file headers.txt
sed -i '1s/^/species\n/' headers.txt
#Combine headers and table with missing data
paste headers.txt missing_percentage_overview.txt > missing_percentage_overview_and_headers.txt
cp missing_percentage_overview_and_headers.txt $path/71selected$CUT

# SELECTION OF MOST COMPLETE ASSEMBLIES
# (i.e., with average percentage of missing data above specified values and with no completely missing data per assembly)
# Transpose data matrix (with percentages of missing data)
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' missing_percentage_overview_and_headers.txt > transposed.txt

# Take first column of space separated file
cut transposed.txt -d' ' -f1 > AssembliesList.txt

# Compute mean from numbers on each line (omitting first column with assembly name), i.e. mean amount of missing data per gene,
# and add this as a last value on each line
awk '
BEGIN {FS=OFS=" "}
{
sum=0; n=0
for(i=1;i<=NF;i++)
     {sum+=$i; ++n}
     print $0," "sum/n
	 sum=0; n=0
}' transposed.txt > transposedPlusMean.txt

# Replace double spaces by single space
sed -i 's/  / /' transposedPlusMean.txt
# Delete first column of a space separated files
awk 'BEGIN{FS=OFS=" "}{$1="";sub(" ","")}1' transposedPlusMean.txt > tmp && mv tmp transposedPlusMean.txt
# Count number of '100' on each line, i.e. number of species with completely missing data, and add this as a last value on each line
awk -F "100" ' { print $0, NF-1 } ' transposedPlusMean.txt > transposedPlusMeanPlusNumberOf100.txt
# Replace two zeros at the end of the first line (result of two previous commands) by column names
# Replace first ' 0 ' on first line with a text ' average_missing_data '
sed -i '1s/ 0 / average_missing_data /' transposedPlusMeanPlusNumberOf100.txt
# Replace first (formely second)' 0' on first line with a text ' nr_assemblies_with_completely_missing_data'
sed -i '1s/ 0/ nr_assemblies_with_completely_missing_data/' transposedPlusMeanPlusNumberOf100.txt
# Combine AssembliesList.txt and transposedPlusMeanPlusNumberOf100.txt
paste AssembliesList.txt transposedPlusMeanPlusNumberOf100.txt > MissingDataOverview.txt
# Print names of assemblies ($1) if mean of missing data (second last column, i.e. $(NF-1)) is lower than specified value
# and number of 100 (i.e., number of completely missing genes, $(NF)) is zero
awk -F' ' -v val=$CUT '( $(NF-1) <= val ) && ( $(NF) < 1 ) { print $1}' MissingDataOverview.txt > selected_genes$CUT.txt

#Copy table and list to home
cp MissingDataOverview.txt $path/71selected${CUT}
cp selected_genes$CUT.txt $path/71selected${CUT}

#Calculate characteristics for selected genes (using AMAS)
mkdir AMAS
cd AMAS
#Copy necessary script
cp $source/AMAS.py .
#Copy selected Assemblies
for file in $(cat ../selected_genes$CUT.txt)
do
	cp ../${file}.fasta .
done
#Calculate summary statistic for selected Assemblies using AMAS
python3 AMAS.py summary -f fasta -d dna -i *.fasta
#Copy summary table to home
cp summary.txt $path/71selected${CUT}

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir04
fi
