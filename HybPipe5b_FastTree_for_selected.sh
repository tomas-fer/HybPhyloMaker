#!/bin/bash
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=4gb
#PBS -N HybPipe5b_FastTree_for_selected
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                   Script 05b - FastTree gene tree building                   *
# *                                                                              *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute gene trees using FastTree for selected genes
# Selection is based on maximum missing data per sample allowed ($MISSINGPERCENT) and minimum species percentage presence per assembly ($SPECIESPRESENCE)
# Edit these values in settings.cfg
# Run first HybPipe4_missingdataremoval.sh with the same settings

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
	module add fasttree-2.1.8
	module add raxml-8.2.4
	module add perl-5.10.1
	module add R-3.1.1
else
	echo "Local..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir05b
	cd workdir05b
fi

#Add necessary scripts and files
cp $source/catfasta2phyml.pl .
cp $source/CompareToBootstrap.pl .
cp $source/MOTree.pm .
cp $path/71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
# Copy and modify selected FASTA files
echo -e "Modifying selested FASTA files...\n"
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
do
	cp $path/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in FastTree)
	sed -i 's/(/_/g' ${i}_modif${MISSINGPERCENT}.fas
	sed -i 's/)//g' ${i}_modif${MISSINGPERCENT}.fas
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	sed -i 's/_contigs//g' ${i}_modif${MISSINGPERCENT}.fas
	sed -i 's/.fas//g' ${i}_modif${MISSINGPERCENT}.fas
done
#Make a list of all fasta files
ls *.fas | cut -d"." -f1 > FileForFastTree.txt
#Make dir for results
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}
mkdir $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree

#----------------Generate gene trees using FastTree----------------
echo -e "Generating FastTrees...\n"
for file in $(cat FileForFastTree.txt)
do
	#FastTree
	if [ ! $LOGNAME == "" ]; then
		echo -e "Processing file: ${file}"
		fasttreemp -nt ${file}.fas > ${file}.fast.tre
	else
		echo -e "Processing file: ${file}"
		fasttree -nt ${file}.fas > ${file}.fast.tre
	fi
	cp *$file.fast.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
done

#----------------Generate bootstrapped gene trees using FastTree----------------
#Bootstrap using FastTree 
echo -e "Generating bootstrapped FastTrees...\n"
if [[ $FastTreeBoot =~ "yes" ]]; then
	for file in $(cat FileForFastTree.txt)
	do
		#Generate 100 replicated datasets using RAxML -f j
		raxmlHPC -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS
		#Loop over replicates and calculate FastTree for each of them
		for i in {0..99}
		do
		if [ ! $LOGNAME == "" ]; then
			fasttreemp -nt ${file}.fas.BS${i} > ${file}.BS${i}.fast.tre
		else
			fasttree -nt ${file}.fas.BS${i} > ${file}.BS${i}.fast.tre
		fi
		done
		#Combine all bootstrap trees to a single file
		cat ${file}.BS*.fast.tre > ${file}.boot.fast.trees
		#Delete BS files and trees
		rm *.BS*
		rm ${file}.BS*.fast.tre
		#Map bootstrap support values onto the original tree
		perl ./CompareToBootstrap.pl -tree *${file}.fast.tre -boot ${file}.boot.fast.trees > ${file}.boot.fast.tre
		#Copy bootstrap trees and a final tree with bootstrap values to home
		cp ${file}.boot.fast.trees $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
		cp ${file}.boot.fast.tre $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	done
fi

#----------------Make a summary table with statistical properties for trees using R----------------
#Copy script
cp $source/tree_props.r .
cp $source/treepropsPlot.r .
cp $source/LBscores.R .
mkdir trees
mkdir alignments
#Copy all fasta alignments to subfolder 'alignments'
cp *.fas alignments/
#Copy all tree files (*.tre) but not bootstrap trees (*boot*) to subfolder 'trees'
cp `find ./*.tre -maxdepth 1 ! -name '*boot*'` trees/
#Run R script for tree properties calculation and for plotting histograms of resulting values
echo -e "\nCalculating tree properties...\n"
R --slave -f tree_props.r
echo -e "Plotting boxplots/histograms for tree properties...\n"
if [ ! $LOGNAME == "" ]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	xvfb-run R --slave -f treepropsPlot.r
else
	R --slave -f treepropsPlot.r
fi
#Run R script for calculation of LB score
echo -e "Calculating and parsing LB score...\n"
R --slave -f LBscores.R
#Parse script output (LBscores.csv)
echo -e "Taxon\tSum\tNrTrees\tMean" > LBscoresPerTaxon.txt
echo -e "Locus\tLBscoreSD" > LBscoresSDPerLocus.txt
for i in $(cat LBscores.csv | sed 1d | cut -d"," -f2 | sort | uniq); do
	awk -F',' -v val=$i 'BEGIN { n=0; sum=0} { if ($2 == val) { sum+=$4; n+=1 } } END { print val "\t" sum "\t" n "\t" sum/n }' LBscores.csv >> LBscoresPerTaxon.txt
done
for i in $(cat LBscores.csv | sed 1d | cut -d"," -f5 | sort | uniq); do
	grep $i LBscores.csv | awk -F',' -v val=$i '{ if ($5 == val) result=$1 } END { print val "\t" result }' >> LBscoresSDPerLocus.txt
done
cp LBscores.csv $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
cp LBscoresPerTaxon.txt $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
cp LBscoresSDPerLocus.txt $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
#Combine 'LBscoresSDPerLocus.txt' with 'tree_stats_table.csv'
awk '{ print $2 }' LBscoresSDPerLocus.txt > tmp && mv tmp LBscoresSDPerLocus.txt
paste tree_stats_table.csv LBscoresSDPerLocus.txt | tr "\t" "," > tmp && mv tmp tree_stats_table.csv
#Copy results to home
cp tree_stats_table.csv $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
cp *.png $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
#Remove results
rm *.png
rm LB*.txt LB*.csv
mv tree_stats_table.csv tree_stats_table_noboot.csv

#----------------Make a summary table for bootstrapped trees----------------
if [[ $FastTreeBoot =~ "yes" ]]; then
	#Remove everything from subfolder 'trees'
	rm ./trees/*.*
	#Copy all trees with bootstrap values to subfolder 'trees'
	cp `find ./*boot.fast.tre -maxdepth 1` trees/
	#Run R script for tree properties calculation
	echo -e "\nCalculating tree properties for bootstrapped trees...\n"
	R --slave -f tree_props.r
	echo -e "Plotting boxplots/histograms for bootstrapped tree properties...\n"
	if [ ! $LOGNAME == "" ]; then
		#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
		xvfb-run R --slave -f treepropsPlot.r
	else
		R --slave -f treepropsPlot.r
	fi
	#Run R script for calculation of LB score
	echo -e "Calculating and parsing LB score...\n"
	R --slave -f LBscores.R
	#Parse script output (LBscores.csv)
	echo -e "Taxon\tSum\tNrTrees\tMean" > LBscoresPerTaxon.txt
	echo -e "Locus\tLBscoreSD" > LBscoresSDPerLocus.txt
	for i in $(cat LBscores.csv | sed 1d | cut -d"," -f2 | sort | uniq); do
		awk -F',' -v val=$i 'BEGIN { n=0; sum=0} { if ($2 == val) { sum+=$4; n+=1 } } END { print val "\t" sum "\t" n "\t" sum/n }' LBscores.csv >> LBscoresPerTaxon.txt
	done
	for i in $(cat LBscores.csv | sed 1d | cut -d"," -f5 | sort | uniq); do
		grep $i LBscores.csv | awk -F',' -v val=$i '{ if ($5 == val) result=$1 } END { print val "\t" result }' >> LBscoresSDPerLocus.txt
	done
	mv LBscores.csv LBscores_bootstrap.csv
	mv LBscoresPerTaxon.txt LBscoresPerTaxon_bootstrap.txt
	mv LBscoresSDPerLocus.txt LBscoresSDPerLocus_bootstrap.txt
	cp LBscores_bootstrap.csv $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	cp LBscoresPerTaxon_bootstrap.txt $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	cp LBscoresSDPerLocus_bootstrap.txt $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	#Combine 'LBscoresSDPerLocus_bootstrap.txt' with 'tree_stats_table.csv'
	awk '{ print $2 }' LBscoresSDPerLocus_bootstrap.txt > tmp && mv tmp LBscoresSDPerLocus_bootstrap.txt
	paste tree_stats_table.csv LBscoresSDPerLocus_bootstrap.txt | tr "\t" "," > tmp && mv tmp tree_stats_table.csv
	#Rename results and copy results to home
	mv tree_stats_table.csv tree_stats_table_bootstrap.csv
	cp tree_stats_table_bootstrap.csv $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	#Rename all PNG files generated by R (add 'bootstrap')and copy them to home
	for file in *.png; do mv "$file" "${file/.png/_bootstrap.png}"; done
	cp *.png $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
fi

#----------------Combine tree summary table with alignment summary and print comparison plots----------------
#Copy script
cp $source/plotting_correlations.R .
echo -e "Combining alignment and tree properties...\n"
#Copy alignment summary
cp $path/71selected${MISSINGPERCENT}/summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
#***Modify alignment summary***
#Remove '_Assembly' (now locus names start with a number)
sed -i 's/Assembly_//g' summarySELECTED*.txt
#Take first line as a header
head -n1 summarySELECTED*.txt > head.txt
#Remove first line | sort
sed 1d summarySELECTED*.txt | sort > summarySELECTED_sorted.txt
#Combine header and sorted file
cat head.txt summarySELECTED_sorted.txt > tmp && mv tmp summarySELECTED_sorted.txt
#***Modify tree summary***
#Change ',' to TAB
sed -i 's/,/\t/g' tree_stats_table_noboot.csv
sed -i 's/,/\t/g' tree_stats_table_bootstrap.csv
#Take first line as a header
head -n1 tree_stats_table_noboot.csv > head2.txt
#Remove first line | sort
sed 1d tree_stats_table_noboot.csv | sort > tree_stats_table_noboot_sorted.csv
sed 1d tree_stats_table_bootstrap.csv | sort > tree_stats_table_bootstrap_sorted.csv
#Combine header and sorted file
cat head2.txt tree_stats_table_noboot_sorted.csv > tmp && mv tmp tree_stats_table_noboot_sorted.csv
cat head2.txt tree_stats_table_bootstrap_sorted.csv > tmp && mv tmp tree_stats_table_bootstrap_sorted.csv
#Combine both files
paste summarySELECTED_sorted.txt tree_stats_table_noboot_sorted.csv > combined_noboot.txt
paste summarySELECTED_sorted.txt tree_stats_table_bootstrap_sorted.csv > combined_bootstrap.txt
#Delete bootstrap-related files if no bootstrap requested
if [[ $FastTreeBoot =~ "no" ]]; then
	rm *bootstrap*
fi
#Rename colums
sed -i 's/Alignment_length/Aln_length/' combined_*.txt
sed -i 's/Missing_percent/Missing_perc/' combined_*.txt
sed -i 's/Proportion_parsimony_informative/Prop_pars_inf/' combined_*.txt
sed -i 's/MstatX_entropy/Aln_entropy/' combined_*.txt
sed -i 's/Average_bootstrap/Bootstrap/' combined_*.txt
sed -i 's/Average_branch_length/Branch_length/' combined_*.txt
sed -i 's/Avg_p_dist/P_distance/' combined_*.txt
sed -i 's/Slope/Satur_slope/' combined_*.txt
sed -i 's/R_squared/Satur_R_sq/' combined_*.txt
sed -i 's/LBscoreSD/LBscore_SD/' combined_*.txt

#Run comparison plots for FastTrees without true bootstrap ('noboot')
mv combined_noboot.txt combined.txt
echo -e "Plotting gene properties correlations for FastTrees without true bootstrap...\n"
if [ ! $LOGNAME == "" ]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	xvfb-run R --slave -f plotting_correlations.R
else
	R --slave -f plotting_correlations.R
fi
cp genes_corrs.* $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
rm genes_corrs.*
mv combined.txt gene_properties.txt
cp gene_properties.txt $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree

#Run comparison plots for FastTrees with true bootstrap ('bootstrap')
echo -e "Plotting gene properties correlations for bootstrapped FastTrees...\n"
if [[ $FastTreeBoot =~ "yes" ]]; then
	mv combined_bootstrap.txt combined.txt
	if [ ! $LOGNAME == "" ]; then
		#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
		xvfb-run R --slave -f plotting_correlations.R
	else
		R --slave -f plotting_correlations.R
	fi
	mv genes_corrs.png genes_corrs_bootstrap.png
	mv genes_corrs.pdf genes_corrs_bootstrap.pdf
	cp genes_corrs_bootstrap.* $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	rm genes_corrs_bootstrap.*
	mv combined.txt gene_properties_bootstrap.txt
	cp gene_properties_bootstrap.txt $path/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
fi

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	#rm -r workdir05b
fi

echo -e "HybPipe 5b finished..."
