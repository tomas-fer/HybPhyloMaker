#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:0:0
#PBS -l select=1:ncpus=4:mem=4gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker6a2_RAxML_trees_summary
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker6a2_RAxML_trees_summary
#$ -o HybPhyloMaker6a2_RAxML_trees_summary.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  Script 06a2 - summary of RAxML gene trees                   *
# *                                   v.1.4.2                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2017 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute summary for already generated gene trees with RAxML
# Run first HybPhyloMaker5_missingdataremoval.sh and HybPhyloMaker6a_RAxML_for_selected.sh with the same settings

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker6a2 is running on MetaCentrum..."
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
	module add R-3.2.3-intel
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker6a2 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir06a2
	cd workdir06a2
	#Add necessary modules
	module load tools/R/3.2.1
else
	echo -e "\nHybPhyloMaker6a2 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir06a2
	cd workdir06a2
fi
#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
else
	echo -e "Working with exons\n"
	type="exons"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
	if [ -d "$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}" ]; then
		if [ "$(ls -A $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT})" ]; then
			if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML" ]; then
				if [ "$(ls -A $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML)" ]; then
					echo -e "OK\n"
				else
					echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' is empty. Exiting...\n"
					rm -d ../workdir06a2/ 2>/dev/null
					exit 3
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' is missing. Exiting...\n"
				rm -d ../workdir06a2/ 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is empty. Exiting...\n"
			rm -d ../workdir06a2/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is missing. Exiting...\n"
		rm -d ../workdir06a2/ 2>/dev/null
		exit 3
	fi
else
	echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
	rm -d ../workdir06a2/ 2>/dev/null
	exit 3
fi

#Test if folder for results exits
if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/gene_properties.txt" ]; then
	echo -e "File '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/gene_properties.txt' already exists. You are probably going to overwrite previous results. Delete this file or rename before running this script again. Exiting...\n"
	rm -d ../workdir06a2/ 2>/dev/null
	exit 3
else
	if [[ ! $location == "1" ]]; then
		if [ "$(ls -A ../workdir06a2)" ]; then
			echo -e "Directory 'workdir06a2' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir06a2/ 2>/dev/null
			exit 3
		fi
	fi
fi

#----------------Make a summary table with statistical properties for trees using R----------------
#Copy script
cp $source/tree_props.r .
cp $source/treepropsPlot.r .
cp $source/LBscores.R .
mkdir trees
mkdir alignments
#Copy all fasta alignments to subfolder 'alignments'
cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt); do
	cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas alignments/
done
#Copy all RAxML tree files (*.tre) with bootstrap values to subfolder 'trees'
cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/*bipartitions.*Assembly* trees/
#Rename RAxML trees
cd trees
for i in *; do
	mv "$i" `basename "$i" .result | cut -f2 -d "."`.tre
done
cd ..

#Set a log file for R outputs/error messages
touch R.log

#Run R script for tree properties calculation and for plotting histograms of resulting values
echo -e "\nCalculating tree properties...\n"
echo -e "Calculating tree properties using tree_props.r\n" >> R.log
R --slave -f tree_props.r >> R.log 2>&1
#Run R script for calculation of LB score
echo -e "Calculating and parsing LB score...\n"
echo -e "\nCalculating and parsing LB score using LBscores.R\n" >> R.log
R --slave -f LBscores.R >> R.log 2>&1
#Parse script output (LBscores.csv)
echo -e "Taxon\tSum\tNrTrees\tMean" > LBscoresPerTaxon.txt
echo -e "Locus\tLBscoreSD" > LBscoresSDPerLocus.txt
for i in $(cat LBscores.csv | sed 1d | cut -d"," -f2 | sort | uniq); do
	awk -F',' -v val=$i 'BEGIN { n=0; sum=0} { if ($2 == val) { sum+=$4; n+=1 } } END { print val "\t" sum "\t" n "\t" sum/n }' LBscores.csv >> LBscoresPerTaxon.txt
done
for i in $(cat LBscores.csv | sed 1d | cut -d"," -f5 | sort | uniq); do
	grep $i LBscores.csv | awk -F',' -v val=$i '{ if ($5 == val) result=$1 } END { print val "\t" result }' >> LBscoresSDPerLocus.txt
done
cp LBscores.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
cp LBscoresPerTaxon.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
cp LBscoresSDPerLocus.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
#Combine 'LBscoresSDPerLocus.txt' with 'tree_stats_table.csv'
awk '{ print $2 }' LBscoresSDPerLocus.txt > tmp && mv tmp LBscoresSDPerLocus.txt
paste tree_stats_table.csv LBscoresSDPerLocus.txt | tr "\t" "," > tmp && mv tmp tree_stats_table.csv
#Replace 'NaN' by '0' (otherwise following plotting in R will not work)
sed -i.bak 's/NaN/0/g' tree_stats_table.csv
echo -e "Plotting boxplots/histograms for tree properties...\n"
echo -e "\nPlotting boxplots/histograms for tree properties using treepropsPlot.r\n" >> R.log
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f treepropsPlot.r
	R --slave -f treepropsPlot.r >> R.log 2>&1
else
	R --slave -f treepropsPlot.r >> R.log 2>&1
fi

#Copy results to home
cp tree_stats_table.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
cp *.png $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

#----------------Combine tree summary table with alignment summary and print comparison plots----------------
#Copy script
cp $source/plotting_correlations.R .
echo -e "Combining alignment and tree properties...\n"
#Copy alignment summary
cp $path/${alnpathselected}${MISSINGPERCENT}/summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
#***Modify alignment summary***
#Remove '_Assembly' (now locus names start with a number)
sed -i.bak 's/Assembly_//g' summarySELECTED*.txt
#Take first line as a header
head -n1 summarySELECTED*.txt > head.txt
#Remove first line | sort
sed 1d summarySELECTED*.txt | sort > summarySELECTED_sorted.txt
#Combine header and sorted file
cat head.txt summarySELECTED_sorted.txt > tmp && mv tmp summarySELECTED_sorted.txt
#***Modify tree summary***
#Change ',' to TAB
sed -i.bak 's/,/\t/g' tree_stats_table.csv
#Take first line as a header
head -n1 tree_stats_table.csv > head2.txt
#Remove first line | sort
sed 1d tree_stats_table.csv | sort > tree_stats_table_sorted.csv
#Combine header and sorted file
cat head2.txt tree_stats_table_sorted.csv > tmp && mv tmp tree_stats_table_sorted.csv
#Combine both files
paste summarySELECTED_sorted.txt tree_stats_table_sorted.csv > combined.txt

#Rename colums
sed -i.bak 's/Alignment_length/Aln_length/' combined.txt
sed -i.bak 's/Missing_percent/Missing_perc/' combined.txt
sed -i.bak 's/Proportion_parsimony_informative/Prop_pars_inf/' combined.txt
sed -i.bak 's/MstatX_entropy/Aln_entropy/' combined.txt
sed -i.bak 's/Average_bootstrap/Bootstrap/' combined.txt
sed -i.bak 's/Average_branch_length/Branch_length/' combined.txt
sed -i.bak 's/Avg_p_dist/P_distance/' combined.txt
sed -i.bak 's/Slope/Satur_slope/' combined.txt
sed -i.bak 's/R_squared/Satur_R_sq/' combined.txt
sed -i.bak 's/LBscoreSD/LBscore_SD/' combined.txt

#Run comparison plots for RAxML trees
echo -e "Plotting gene properties correlations for RAxML trees...\n"
echo -e "\nPlotting gene properties correlations for RAxML trees using plotting_correlations.R\n" >> R.log
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f plotting_correlations.R
	R --slave -f plotting_correlations.R >> R.log 2>&1
else
	R --slave -f plotting_correlations.R >> R.log 2>&1
fi
cp genes_corrs.* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
rm genes_corrs.*
mv combined.txt gene_properties.txt
cp gene_properties.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

#Copy R.log to home
cp R.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir06a2
fi

echo -e "HybPhyloMaker 6a2 finished...\n"
