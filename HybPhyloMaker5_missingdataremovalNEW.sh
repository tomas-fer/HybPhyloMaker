#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=2gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker5_missingdata_handling
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q mThC.q
#$ -l mres=2G,h_data=2G,h_vmem=2G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker5_missingdata_handling
#$ -o HybPhyloMaker5_missingdata_handling.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                      Script 05 - Missing data handling                       *
# *                                   v.1.6.7                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2018 *
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
	echo -e "\nHybPhyloMaker5 is running on MetaCentrum..."
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
	module add seqtk-1.0
	module add R-3.4.3-gcc
	module add debian8-compat
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker5 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir05
	cd workdir05
	#Add necessary modules
	module load bioinformatics/anaconda3/2.3.0 #for python3
	module load bioinformatics/trimal/1.4
	module load bioinformatics/mstatx/1.0
	module load tools/R/3.4.1
else
	echo -e "\nHybPhyloMaker5 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	#echo -e "Path is: $path"
	source=../HybSeqSource
	#echo -e "Source is: $source\n"
	#Make and enter work directory
	mkdir -p workdir05
	cd workdir05
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
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
	echo -e "...with corrected reading frame\n"
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
	echo -e "\n"
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [[ $cp =~ "yes" ]]; then
	if [ -d "$path/$type/60mafft" ]; then
		if [ "$(ls -A $path/$type/60mafft)" ]; then
			echo -e "OK\n"
		else
			echo -e "'$path/$type/60mafft' is empty. Exiting...\n"
			rm -d ../workdir05/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/$type/60mafft' is missing. Exiting...\n"
		rm -d ../workdir05/ 2>/dev/null
		exit 3
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		if [ -d "$path/$type/80concatenated_exon_alignments_corrected" ]; then
			if [ "$(ls -A $path/$type/80concatenated_exon_alignments_corrected)" ]; then
				echo -e "OK\n"
			else
				echo -e "'$path/$type/80concatenated_exon_alignments_corrected' is empty. Exiting...\n"
				rm -d ../workdir05/ 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/$type/80concatenated_exon_alignments_corrected' is missing. Exiting...\n"
			rm -d ../workdir05/ 2>/dev/null
			exit 3
		fi
	else
		if [ -d "$path/$type/70concatenated_exon_alignments" ]; then
			if [ "$(ls -A $path/$type/70concatenated_exon_alignments)" ]; then
				echo -e "OK\n"
			else
				echo -e "'$path/$type/70concatenated_exon_alignments' is empty. Exiting...\n"
				rm -d ../workdir05/ 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/$type/70concatenated_exon_alignments' is missing. Exiting...\n"
			rm -d ../workdir05/ 2>/dev/null
			exit 3
		fi
	fi
fi

#Test if folder for results exits
if [[ $corrected =~ "yes" ]]; then
	if [ -d "$path/$type/81selected_corrected${MISSINGPERCENT}" ]; then
		echo -e "Directory '$path/$type/81selected_corrected${MISSINGPERCENT}' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir05/ 2>/dev/null
		exit 3
	else
		if [[ ! $location == "1" ]]; then
			if [ "$(ls -A ../workdir05)" ]; then
				echo -e "Directory 'workdir05' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir05/ 2>/dev/null
				exit 3
			fi
		fi
	fi
else
	if [ -d "$path/$type/71selected${MISSINGPERCENT}" ]; then
		echo -e "Directory '$path/$type/71selected${MISSINGPERCENT}' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir05/ 2>/dev/null
		exit 3
	else
		if [[ ! $location == "1" ]]; then
			if [ "$(ls -A ../workdir05)" ]; then
				echo -e "Directory 'workdir05' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir05/ 2>/dev/null
				exit 3
			fi
		fi
	fi
fi
#Copy data folder to scratch
if [[ $cp =~ "yes" ]]; then
	cp $path/$type/60mafft/*.fasta .
	#Rename *.mafft to *.fasta
	# for file in *.mafft; do
		# mv "$file" "${file%.fasta.mafft}.fasta"
	# done
else
	#cp $path/$alnpath/*.fasta .
	find $path/$alnpath -maxdepth 1 -name "*.fasta" -exec cp -t . {} + #to avoid 'Argument list too long' error
fi

#-----------PREPARE ALIGNMENTS WITH SPECIES WITH MAXIMUM SPECIFIED MISSING DATA ONLY----------------------
#Make a list of all fasta files
#ls *.fasta | cut -d"." -f1 > fileForDeletePercentage.txt
find -name "*.fasta" | sed 's/\.\///' | sed 's/\.fasta//g' > fileForDeletePercentage.txt #to avoid 'Argument list too long' error
#Make new dir for results
mkdir $path/${alnpathselected}${MISSINGPERCENT}
mkdir $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
numberfiles=$(cat fileForDeletePercentage.txt | wc -l)
calculating=0
for file in $(cat fileForDeletePercentage.txt)
do
	#Delete empty lines (to be on the safe side...)
	sed -i.bak '/^$/d' $file.fasta
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file.fasta > tmp && mv tmp $file.fasta
	#modify fasta
	#1. replace leading and tailing '-' to '?'
	perl -pe 's/\G-|-(?=-*$)/?/g' $file.fasta > ${file}_modif.fasta
	#2. change all valid characters to 'X'
	sed -i.bak '/^>/!s/[ACGTRYSWKMBDHVacgtryswkmbdhv]/X/g' ${file}_modif.fasta
	#3. change missing to 'A'
	sed -i.bak2 '/^>/!s/[nN?]/A/g' ${file}_modif.fasta
	#make table with length, number and percentage of missing data
	seqtk comp ${file}_modif.fasta | awk '{print $1,$3/$2*100}' > ${file}_${MISSINGPERCENT}percN.fas
	#make a list of sequences to keep (i.e., with missing data below the threshold)
	awk -v val=$MISSINGPERCENT '$2<=val {print $1}' ${file}_${MISSINGPERCENT}percN.fas > keep.txt
	#add header to the file with missing percent
	sed -i.bak "1s/^/species $file\n/" ${file}_${MISSINGPERCENT}percN.fas
	#change separator from ' ' to TABs
	tr ' ' '\t' < ${file}_${MISSINGPERCENT}percN.fas > tmp && mv tmp ${file}_${MISSINGPERCENT}percN.fas
	#grep only sequences from keep.txt (and remove separators and empty lines)
	grep -f keep.txt -A1 $file.fasta | sed 's/^--$//' | sed '/^$/d'> ${file}_modif${MISSINGPERCENT}.fas
	#Copy results home
	cp ${file}_${MISSINGPERCENT}percN.fas $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
	cp ${file}_modif${MISSINGPERCENT}.fas $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
	calculating=$((calculating + 1))
	rm ${file}_modif.fasta keep.txt
	echo -e "\t$file processed ($calculating out of $numberfiles)"
done

echo -e "\nDeleting samples with more than ${MISSINGPERCENT}% missing data finished."
echo -e "Modified alignments saved in $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}"

#-----------MAKE A TABLE WITH % OF MISSING DATA IN EACH SPECIES AND ASSEMBLY----------------------
#Prepare header file
echo -e "\nPreparing summary table..."
file=$(cat fileForDeletePercentage.txt | head -n 1)
cat $file.fasta | sed '/^>/{N; s/\n/ /;}' | cut -f1 -d" " | sed 's/>//' | sed 's/_contigs//' | sed 's/.fas//' > headers.txt

#Make a table with percentages of N in each accession and file
#awk '{_[FNR]=(_[FNR] OFS $2)}END{for (i=1; i<=FNR; i++) {sub(/^ /,"",_[i]); print _[i]}}' *percN.fas > missing_percentage_overview.txt
find -name "*percN.fas" -print0 | xargs -0 awk '{_[FNR]=(_[FNR] OFS $2)}END{for (i=1; i<=FNR; i++) {sub(/^ /,"",_[i]); print _[i]}}' > missing_percentage_overview.txt #to avoid 'Argument list too long' error
#Add word 'species' as a 1st line in file headers.txt
sed -i.bak '1s/^/species\n/' headers.txt
#Combine headers and table with missing data
paste headers.txt missing_percentage_overview.txt > missing_percentage_overview_and_headers.txt
cp missing_percentage_overview_and_headers.txt $path/${alnpathselected}${MISSINGPERCENT}

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
sed -i.bak 's/  / /' transposedPlusMean.txt
# Delete first column of a space separated files
awk 'BEGIN{FS=OFS=" "}{$1="";sub(" ","")}1' transposedPlusMean.txt > tmp && mv tmp transposedPlusMean.txt
# Count number of '100' on each line, i.e. number of species with completely missing data, and add this as a last value on each line
awk -F "100" ' { print $0, NF-1 } ' transposedPlusMean.txt > transposedPlusMeanPlusNumberOf100.txt
# Replace two zeros at the end of the first line (result of two previous commands) by column names
# Replace first ' 0 ' on first line with a text ' average_missing_data '
sed -i.bak '1s/ 0 / average_missing_data /' transposedPlusMeanPlusNumberOf100.txt
# Replace first (formely second)' 0' on first line with a text ' nr_assemblies_with_completely_missing_data'
sed -i '1s/ 0/ nr_assemblies_with_completely_missing_data/' transposedPlusMeanPlusNumberOf100.txt
# Combine AssembliesList.txt and transposedPlusMeanPlusNumberOf100.txt
paste AssembliesList.txt transposedPlusMeanPlusNumberOf100.txt > MissingDataOverview.txt
# Copy table and list to home
cp MissingDataOverview.txt $path/${alnpathselected}${MISSINGPERCENT}
echo -e "Table with % of missing data per gene and sample saved to $path/$type/71selected${MISSINGPERCENT}\n"

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
sed -i.bak '1d' MissingDataOverview_${MISSINGPERCENT}.txt
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
# Select assemblies with equal or more than SPECIESPRESENCE of species (and delete first and last line including header and sum legend)
awk -F' ' -v val=$SPECIESPRESENCE '( $(NF) >= val/100 ) { print $1 }' MissingDataOverview_${MISSINGPERCENT}.txt | tail -n +2 | head -n -2 > selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
# Copy table and list to home
cp MissingDataOverview_${MISSINGPERCENT}.txt $path/${alnpathselected}${MISSINGPERCENT}
cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt $path/${alnpathselected}${MISSINGPERCENT}
echo -e "Assemblies including at least ${SPECIESPRESENCE}% of species selected"

#-----------CALCULATE CHARACTERISTICS FOR ALL GENES AND FOR SELECTED GENES (USING AMAS)----------------------
# Copy necessary script
cp $source/AMAS.py .
cp $source/alignmentSummary.R .
# 1. For all genes
#Calculate alignment summary using AMAS
echo -e "\nCalculating alignment characteristics for all genes using AMAS..."
#This is a faster solution but with really many genes generate 'Argument list too long' error
#python3 AMAS.py summary -f fasta -d dna -i *.fasta

#Much slower solution but works also in case of many genes
for f in *.fasta ; do python3 AMAS.py summary -f fasta -d dna -i $f -o ${f}.summary; done #one summary per gene
find -name "*.summary" -print0 | xargs -0 awk 'FNR==2' > amas.sum #combine all *.summary (take 2nd line, i.e. the data, from each summary by AMAS)
head -n1 `find -name "*.summary" -print -quit` > amas.header #make header (take first line from the first match to *.summary)
cat amas.header amas.sum > summary.txt

#Calculate global alignment entropy using MstatX
echo -e "\nCalculating alignment entropy for all genes using MstatX..."
for file in $(ls *.fasta); do
	echo $file
	mstatx -i $file -g > /dev/null 2>&1
	line=$(cat output.txt)
	echo -e "$file\t$line" >> mstatx.txt
done
rm output.txt
#Add header to the first line
#sed -i.bak -e "1iLocus\tMstatX_entropy" mstatx.txt
echo -e "Locus\tMstatX_entropy" > mstatx.header
cat mstatx.header mstatx.txt > tmp && mv tmp mstatx.txt
#Take only 2nd column with MstatX results (omitting alignment name)
awk '{ print $2 }' mstatx.txt > tmp && mv tmp mstatx.txt
#Replace '?' by 'n' and all 'n' by 'N'
#sed -i.bak 's/\?/n/g' *.fasta
find -name "*fasta" -print0 | xargs -0 sed -i.bak 's/\?/n/g' #to avoid 'Argument list too long' error
#sed -i.bak '/^>/!s/n/N/g' *.fasta
find -name "*fasta" -print0 | xargs -0 sed -i.bak '/^>/!s/n/N/g' #to avoid 'Argument list too long' error
#Calculate alignment conservation value using trimAl
echo -e "\nCalculating alignment conservation value for all genes using trimAl..."
for file in $(ls *.fasta); do
	echo $file
	#Calculate sct using trimal | delete first three lines (header) | replace double TABs by single | replace ' ' by nothing | calculate weighted mean by
	#multiplying first column value (number of residues, i.e. positions) by fifth column (similarity value), sum up all values and divide by the length of the gene
	trimal -in $file -sst | sed '1,3d' | sed 's/\t\t/\t/g' | sed 's/ //g' | awk -v val=$file '{ sum+=$1*$5; sum2+=$1} END { print val "\t" sum/sum2}' >> sct.out
done
#Add header to the first line
#sed -i.bak -e "1iLocus\ttrimAl_sct" sct.out
echo -e "Locus\ttrimAl_sct" > trimal.header
cat trimal.header sct.out > tmp && mv tmp sct.out
#Take only 2nd column with trimAl results (omitting alignment name)
awk '{ print $2 }' sct.out > tmp && mv tmp sct.out
#Combine results together
paste summary.txt mstatx.txt sct.out > summaryALL.txt
#Correct '-inf' (generated by TrimAl when no variation) to '1'
sed -i 's/-inf/1/g' summaryALL.txt
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
cp summaryALL.txt $path/${alnpathselected}${MISSINGPERCENT}
cp *.png $path/${alnpathselected}${MISSINGPERCENT}
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
#This is a faster solution but with really many genes generate 'Argument list too long' error
#python3 AMAS.py summary -f fasta -d dna -i AMASselected/*.fas

#Much slower solution but works also in case of many genes
for f in AMASselected/*.fas ; do python3 AMAS.py summary -f fasta -d dna -i $f -o ${f}.summary; done #one summary per gene
find AMASselected/ -name "*.summary" -print0 | xargs -0 awk 'FNR==2' > amasSel.sum #combine all *.summary (take 2nd line, i.e. the data, from each summary by AMAS)
head -n1 `find AMASselected/ -name "*.summary" -print -quit` > amasSel.header #make header (take first line from the first match to *.summary)
cat amasSel.header amasSel.sum > summary.txt

#Calculate global alignment entropy using MstatX
echo -e "\nCalculating alignment entropy for selected genes using MstatX..."
for file in $(ls AMASselected/*.fas); do
	echo $file
	mstatx -i $file -g  > /dev/null 2>&1
	line=$(cat output.txt)
	echo -e "$file\t$line" >> mstatx.txt
done
rm output.txt
#Add header to the first line
sed -i.bak -e "1iLocus\tMstatX_entropy" mstatx.txt
#Take only 2nd column with MstatX results (omitting alignment name)
awk '{ print $2 }' mstatx.txt > tmp && mv tmp mstatx.txt
#Replace '?' by 'n'
#sed -i.bak 's/\?/n/g' AMASselected/*.fas
find AMASselected/ -name "*fas" -print0 | xargs -0 sed -i.bak 's/\?/n/g' #to avoid 'Argument list too long' error
#sed -i.bak 's/n/N/g' AMASselected/*.fas
find AMASselected/ -name "*fas" -print0 | xargs -0 sed -i.bak '/^>/!s/n/N/g' #to avoid 'Argument list too long' error
#Calculate conservation value using trimAl
echo -e "\nCalculating alignment conservation value for selected genes using trimAl..."
for file in $(ls AMASselected/*.fas); do
	echo $file
	#Calculate sct using trimal | delete first three lines (header) | replace double TABs by single | replace ' ' by nothing | calculate weighted mean by
	#multiplying first column value (number of residues, i.e. positions) by fifth column (similarity value), sum up all values and divide by the length of the gene
	trimal -in $file -sst | sed '1,3d' | sed 's/\t\t/\t/g' | sed 's/ //g' | awk -v val=$file '{ sum+=$1*$5; sum2+=$1} END { print val "\t" sum/sum2}' >> sct.out
done
#Add header to the first line
sed -i.bak -e "1iLocus\ttrimAl_sct" sct.out
#Take only 2nd column with trimAl results (omitting alignment name)
awk '{ print $2 }' sct.out > tmp && mv tmp sct.out
#Combine results together
paste summary.txt mstatx.txt sct.out > summaryALL.txt
#Correct '-inf' (generated by TrimAl when no variation) to '1'
sed -i 's/-inf/1/g' summaryALL.txt
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
cp summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt $path/${alnpathselected}${MISSINGPERCENT}
cp *.png $path/${alnpathselected}${MISSINGPERCENT}

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir05
fi

echo -e "\nHybPhyloMaker5 finished...\n"
