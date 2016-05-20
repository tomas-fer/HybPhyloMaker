#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=12gb
#PBS -l scratch=4gb
#PBS -N HybPipe4_missingdata_handling
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPipe4_missingdata_handling
#$ -o HybPipe4_missingdata_handling.log


# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                      Script 04 - Missing data handling                       *
# *                                   v.1.0.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Remove species with given percentage (or more) of missing data specified in $MISSINGPERCENT (i.e., 'N', 'n' or '?', but NOT '-') from multiple FASTA file
#Removal is done for all *.fasta files in the folder 70concatenated_exon_alignments
#All fasta files must include all samples in the same order!!!
#Headers must not include spaces!!!
#Creates following files:
#(1) modified FASTA files without samples with percentage of missing data higher than specified value
#(2) table with headers and percentage of missing data before removal of samples with excessive amount of missing data
#(3) table with headers and percentage of missing data after removal
#(4) list of assemblies with data for at least $SPECIESPRESENCE percentage of samples
#(5) table of summary characteristics (length, %missing, number of variable and parsimony informative sites, GC content...) for all loci before removal and each selected locus after removal
#(6) histograms/boxplots of selected summary characteristic for all loci before removal and each selected locus after removal (using external R script)

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
	#Add necessary modules
	module add python-3.4.1-gcc
	module add trimal-1.4
	module add mstatx
	module add R-3.2.3-intel
elif [[ $HOSTNAME == *local* ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir04
	cd workdir04
	#Add necessary modules
	module load bioinformatics/anaconda3/2.3.0
	module load tools/R/3.2.1
	
else
	echo "Local..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	echo -e "Path is: $path"
	source=../HybSeqSource
	echo -e "Source is: $source\n"
	#Make and enter work directory
	mkdir workdir04
	cd workdir04
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	type="_cp"
else
	type=""
fi

#Copy data folder to scratch
if [[ $cp =~ "yes" ]]; then
	cp $path/60mafft_cp/*.fasta .
	#Rename *.mafft to *.fasta
	# for file in *.mafft; do
		# mv "$file" "${file%.fasta.mafft}.fasta"
	# done
else
	cp $path/70concatenated_exon_alignments/*.fasta .
fi

#-----------PREPARE ALIGNMENTS WITH SPECIES WITH MAXIMUM SPECIFIED MISSING DATA ONLY----------------------
#Make a list of all fasta files
ls *.fasta | cut -d"." -f1 > fileForDeletePercentage.txt
#Make new dir for results
mkdir $path/71selected${type}${MISSINGPERCENT}
mkdir $path/71selected${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
for file in $(cat fileForDeletePercentage.txt)
do
	#Delete empty lines (to be on the safe side...)
	sed -i '/^$/d' $file.fasta
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file.fasta > tmp && mv tmp $file.fasta
	#Replace newline with ' ' if line starts with '>' (i.e., merge headers with data into single line separated by space)
	sed -i '/^>/{N; s/\n/ /;}' $file.fasta
	
	#Create file for storage headers and percentages of N
	touch ${file}_${MISSINGPERCENT}percN.fas
	echo -e "species\t$file" >> ${file}_${MISSINGPERCENT}percN.fas
	#Create file for accessions with N below threshold (retained accessions)
	touch ${file}_modif${MISSINGPERCENT}.fas
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
		echo -e "$header\t$percentageprint" >> ${file}_${MISSINGPERCENT}percN.fas
		if [ "$percentage" -lt "$MISSINGPERCENT" ]; then
			echo "$line" >> ${file}_modif${MISSINGPERCENT}.fas
		else
			echo -e "Deleting $header with $percentageprint of missing data"
		fi
	done
	#Replace ' ' by newline character
	sed -i '/^>/{s/ /\n/}' $file.fasta
	sed -i '/^>/{s/ /\n/}' ${file}_modif${MISSINGPERCENT}.fas
	#Copy results home
	cp ${file}_${MISSINGPERCENT}percN.fas $path/71selected${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
	cp ${file}_modif${MISSINGPERCENT}.fas $path/71selected${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
	echo -e "\t$file processed\n"
done
echo -e "\nDeleting samples with more than ${MISSINGPERCENT}% missing data finished."
echo -e "Modified alignments saved in $path/71selected${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}"

#-----------MAKE A TABLE WITH % OF MISSING DATA IN EACH SPECIES AND ASSEMBLY----------------------
#Prepare header file
file=$(cat fileForDeletePercentage.txt | head -n 1)
cat $file.fasta | sed '/^>/{N; s/\n/ /;}' | cut -f1 -d" " | sed 's/>//' | sed 's/_contigs//' | sed 's/.fas//' > headers.txt

#Make a table with percentages of N in each accession and file
awk '{_[FNR]=(_[FNR] OFS $2)}END{for (i=1; i<=FNR; i++) {sub(/^ /,"",_[i]); print _[i]}}' *percN.fas > missing_percentage_overview.txt
#Add word 'species' as a 1st line in file headers.txt
sed -i '1s/^/species\n/' headers.txt
#Combine headers and table with missing data
paste headers.txt missing_percentage_overview.txt > missing_percentage_overview_and_headers.txt
cp missing_percentage_overview_and_headers.txt $path/71selected${type}${MISSINGPERCENT}

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
# Copy table and list to home
cp MissingDataOverview.txt $path/71selected${type}${MISSINGPERCENT}
echo -e "Table with % of missing data per gene and sample saved to $path/71selected${MISSINGPERCENT}\n"

#-----------SELECTION OF MOST COMPLETE ASSEMBLIES----------------------
# (i.e., only containing species with at least $MISSINGPERCENT data and are present at least in $SPECIESPRESENCE of samples)

# This selection is no longer used..., SEE FURTHER
# i.e., WE DO NOT USE selection of genes with average percentage of missing data above specified values and with no completely missing data per assembly
# Print names of assemblies ($1) if mean of missing data (second last column, i.e. $(NF-1)) is lower than specified value
# and number of 100 (i.e., number of completely missing genes, $(NF)) is zero
# awk -F' ' -v val=$CUT '( $(NF-1) <= val ) && ( $(NF) < 1 ) { print $1}' MissingDataOverview.txt > selected_genes$CUT.txt
# cp selected_genes$CUT.txt $path/71selected${MISSINGPERCENT}

echo "Selection of most complete assemblies..."
# Remove last two columns from MissingDataOverview.txt
cat MissingDataOverview.txt | awk 'NF{NF-=2}1' FS=' ' OFS=' ' > MissingDataOverview_${MISSINGPERCENT}.txt
# Copy first line to header.txt
cat MissingDataOverview_${MISSINGPERCENT}.txt | head -n1 > header.txt
# Remove first line
sed -i '1d' MissingDataOverview_${MISSINGPERCENT}.txt
# Count nr. of loci (=i.e., nr. of lines)
nrgenes=$(cat MissingDataOverview_${MISSINGPERCENT}.txt | wc -l)
echo "Number of genes: " $nrgenes
# Replace all values higher than MISSINGPERCENT by N/A (except in the first column, i.e. starting in 2nd column)
awk -F" " -v val=$MISSINGPERCENT '{for(i=2;i<=NF;i++) if($i>val) {$i="N/A"}}1' OFS=" " MissingDataOverview_${MISSINGPERCENT}.txt > tmp && mv tmp MissingDataOverview_${MISSINGPERCENT}.txt
# Count number of species (columns - 1)
nrspecies=$(awk '{print NF-1; exit}' MissingDataOverview_${MISSINGPERCENT}.txt)
echo "Number of species: " $nrspecies
# Compute average from numbers on each line (omitting first column with assembly name), i.e. average amount of missing data per gene (after N/A replacement)
# and add this as a last value on each line
awk '{sum=cnt=0; for (i=2;i<=NF;i++) if ($i != "N/A") { sum+=$i; cnt++ } print $0, (cnt ? sum/cnt : "N/A") }' MissingDataOverview_${MISSINGPERCENT}.txt > tmp && mv tmp MissingDataOverview_${MISSINGPERCENT}.txt
# Count number of 'N/A' on each line, i.e. number of species with completely missing data, calculate percentage of non-missing entries and add this as a last value on each line
awk -F "N/A" -v val=$nrspecies ' { print $0, (val-(NF-1))/val } ' MissingDataOverview_${MISSINGPERCENT}.txt > tmp && mv tmp MissingDataOverview_${MISSINGPERCENT}.txt
# Add two last column to header.txt
awk '{ print $0 " averageMissing percPresentSpecies"}' header.txt > tmp && mv tmp header.txt
# Combine header and missing overview
cat header.txt MissingDataOverview_${MISSINGPERCENT}.txt > tmp && mv tmp MissingDataOverview_${MISSINGPERCENT}.txt
# Count number of 'N/A' and calculate average in each column
for i in $(seq 2 `expr $nrspecies + 1`); do
	# Number of genes with 'N/A'
	geneabsence=$(cat MissingDataOverview_${MISSINGPERCENT}.txt | cut -d" " -f$i | grep N/A | wc -l)
	# Number of genes without 'N/A' (i.e., number of genes with at least $MISSINGPERCENT data
	genepresence=`expr $nrgenes - $geneabsence`
	echo $genepresence >> genespresent_${MISSINGPERCENT}.txt
	# Take i-th column, delete two last lines, remove N/A, delete empty lines, calculate average from remaining values
	cat MissingDataOverview_${MISSINGPERCENT}.txt | cut -d" " -f$i | tail -n +2 | sed 's/N\/A//g' | sed '/^\s*$/d' | awk '{ sum += $1 } END { print sum / NR }' >> average_missing_${MISSINGPERCENT}.txt
done
# Make last two lines (with total counts of 'N/A' and with average missing data)
cat genespresent_${MISSINGPERCENT}.txt | tr "\n" " " | awk '{ print "total_genes " $0 }' > tmp && mv tmp genespresent_${MISSINGPERCENT}.txt
cat average_missing_${MISSINGPERCENT}.txt | tr "\n" " " | awk '{ print "average_missing " $0 }' > tmp && mv tmp average_missing_${MISSINGPERCENT}.txt
cat MissingDataOverview_${MISSINGPERCENT}.txt average_missing_${MISSINGPERCENT}.txt genespresent_${MISSINGPERCENT}.txt > tmp && mv tmp MissingDataOverview_${MISSINGPERCENT}.txt
# Select assemblies with more than 75% of species (and delete first and last line including header and sum legend)
awk -F' ' -v val=$SPECIESPRESENCE '( $(NF) > val/100 ) { print $1 }' MissingDataOverview_${MISSINGPERCENT}.txt | tail -n +2 | head -n -2 > selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
# Copy table and list to home
cp MissingDataOverview_${MISSINGPERCENT}.txt $path/71selected${type}${MISSINGPERCENT}
cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt $path/71selected${type}${MISSINGPERCENT}
echo -e "Assemblies including at least ${SPECIESPRESENCE}% of  species selected"

#-----------CALCULATE CHARACTERISTICS FOR ALL GENES AND FOR SELECTED GENES (USING AMAS)----------------------
# Copy necessary script
cp $source/AMAS.py .
cp $source/alignmentSummary.R .
# 1. For all genes
#Calculate alignment summary using AMAS
echo -e "\nCalculating alignment characteristics for all genes using AMAS..."
python3 AMAS.py summary -f fasta -d dna -i *.fasta
#Calculate global alignment entropy using MstatX
echo -e "\nCalculating alignment entropy for all genes using MstatX..."
for file in $(ls *.fasta); do
	mstatx -i $file -g
	line=$(cat output.txt)
	echo -e "$file\t$line" >> mstatx.txt
done
rm output.txt
#Add header to the first line
sed -i -e "1iLocus\tMstatX_entropy" mstatx.txt
#Take only 2nd column with MstatX results (omitting alignment name)
awk '{ print $2 }' mstatx.txt > tmp && mv tmp mstatx.txt
#Replace '?' by 'n' and all 'n' by 'N'
sed -i 's/\?/n/g' *.fasta
sed -i 's/n/N/g' *.fasta
#Calculate alignment conservation value using trimAl
echo -e "Calculating alignment conservation value for all genes using trimAl..."
for file in $(ls *.fasta); do
	echo $file
	#Calculate sct using trimal | delete first three lines (header) | replace double TABs by single | replace ' ' by nothing | calculate weighted mean by
	#multiplying first column value (number of residues, i.e. positions) by fifth column (similarity value), sum up all values and divide by the length of the gene
	trimal -in $file -sst | sed '1,3d' | sed 's/\t\t/\t/g' | sed 's/ //g' | awk -v val=$file '{ sum+=$1*$5; sum2+=$1} END { print val "\t" sum/sum2}' >> sct.out
done
#Add header to the first line
sed -i -e "1iLocus\ttrimAl_sct" sct.out
#Take only 2nd column with trimAl results (omitting alignment name)
awk '{ print $2 }' sct.out > tmp && mv tmp sct.out
#Combine results together
paste summary.txt mstatx.txt sct.out > summaryALL.txt
rm summary.txt mstatx.txt sct.out
echo -e "\nPlotting boxplots/histograms for alignment characteristics for all genes ..."
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f alignmentSummary.R
	R --slave -f alignmentSummary.R
else
	R --slave -f alignmentSummary.R
fi
#Copy summary table to home
cp summaryALL.txt $path/71selected${type}${MISSINGPERCENT}
cp *.png $path/71selected${type}${MISSINGPERCENT}
rm summaryALL.txt
rm *.png
# 2. For selected genes
mkdir AMASselected
#Copy selected Assemblies
for file in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
do
	cp ${file}_modif${MISSINGPERCENT}.fas AMASselected/
done
#Calculate summary statistic for selected Assemblies using AMAS
echo -e "\nCalculating alignment characteristics for selected genes using AMAS..."
python3 AMAS.py summary -f fasta -d dna -i AMASselected/*.fas
#Calculate global alignment entropy using MstatX
echo -e "\nCalculating alignment entropy for selected genes using MstatX..."
for file in $(ls AMASselected/*.fas); do
	mstatx -i $file -g
	line=$(cat output.txt)
	echo -e "$file\t$line" >> mstatx.txt
done
rm output.txt
#Add header to the first line
sed -i -e "1iLocus\tMstatX_entropy" mstatx.txt
#Take only 2nd column with MstatX results (omitting alignment name)
awk '{ print $2 }' mstatx.txt > tmp && mv tmp mstatx.txt
#Replace '?' by 'n'
sed -i 's/\?/n/g' AMASselected/*.fas
sed -i 's/n/N/g' AMASselected/*.fas
#Calculate conservation value using trimAl
echo -e "Calculating alignment conservation value for selected genes using trimAl..."
for file in $(ls AMASselected/*.fas); do
	echo $file
	#Calculate sct using trimal | delete first three lines (header) | replace double TABs by single | replace ' ' by nothing | calculate weighted mean by
	#multiplying first column value (number of residues, i.e. positions) by fifth column (similarity value), sum up all values and divide by the length of the gene
	trimal -in $file -sst | sed '1,3d' | sed 's/\t\t/\t/g' | sed 's/ //g' | awk -v val=$file '{ sum+=$1*$5; sum2+=$1} END { print val "\t" sum/sum2}' >> sct.out
done
#Add header to the first line
sed -i -e "1iLocus\ttrimAl_sct" sct.out
#Take only 2nd column with trimAl results (omitting alignment name)
awk '{ print $2 }' sct.out > tmp && mv tmp sct.out
#Combine results together
paste summary.txt mstatx.txt sct.out > summaryALL.txt
rm summary.txt mstatx.txt sct.out
echo -e "\nPlotting boxplots/histograms for alignment characteristics for selected genes ..."
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	xvfb-run R --slave -f alignmentSummary.R
else
	R --slave -f alignmentSummary.R
fi
mv summaryALL.txt summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#Rename all PNG files generated by R (add '_${MISSINGPERCENT}_${SPECIESPRESENCE}')
for file in *.png; do mv "$file" "${file/.png/_${MISSINGPERCENT}_${SPECIESPRESENCE}.png}"; done
#Copy summary table and PNG files to home
cp summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt $path/71selected${type}${MISSINGPERCENT}
cp *.png $path/71selected${type}${MISSINGPERCENT}
echo -e "\nHybPipe4 finished..."

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir04
fi
