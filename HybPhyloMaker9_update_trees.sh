#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker9_update_trees
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker9_update_trees
#$ -o HybPhyloMaker9_update_trees.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                           Script 09 - Update trees                           *
# *                                   v.1.5.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2017 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#UPDATE on tree selection
# 

if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker9 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add R-3.2.3-intel
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker9 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir09
	cd workdir09
	#Add necessary modules
	module load tools/R/3.2.1
else
	echo -e "\nHybPhyloMaker9 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir09
	cd workdir09
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
	echo -en "...with corrected reading frame"
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
fi

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection\n"
else
	echo -e "\n"
fi

#Copy scripts
cp $source/plotting_correlations.R .
cp $source/alignmentSummary.R .
cp $source/treepropsPlot.r .
#Set a log file for R outputs/error messages
touch R.log

#Copy updated gene list with properties (with unwanted genes deleted)
cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/gene_properties_update.txt .
#Test if there is enough genes selected (i.e., 10 or more)
if [[ $(cat gene_properties_update.txt | wc -l) -lt 11 ]]; then
	echo -e "Less than 10 genes selected. Exiting...\n"
	cd ..
	rm -r workdir09
	exit 3
fi

#Plot gene properties correlations
cp gene_properties_update.txt combined.txt
echo -e "Plotting gene properties correlations for updated selection..."
echo -e "\nPlotting gene properties correlations for updated selection using plotting_correlations.R\n" >> R.log
if [ ! $location == "1" ]; then
	#Run R script for correlation visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#echo "Running xvfb-run R..."
	#xvfb-run R --slave -f plotting_correlations.R
	R --slave -f plotting_correlations.R >> R.log 2>&1
else
	R --slave -f plotting_correlations.R >> R.log 2>&1
fi
mv genes_corrs.png genes_corrs_update.png
mv genes_corrs.pdf genes_corrs_update.pdf
cp genes_corrs_update.* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/
#Plot boxplots/histograms for selected alignment properties
cp gene_properties_update.txt summaryALL.txt
sed -i.bak 's/Aln_length/Alignment_length/' summaryALL.txt
sed -i.bak 's/Missing_perc/Missing_percent/' summaryALL.txt
sed -i.bak 's/Prop_pars_inf/Proportion_parsimony_informative/' summaryALL.txt
sed -i.bak 's/Aln_entropy/MstatX_entropy/' summaryALL.txt

echo -e "\nPlotting boxplots/histograms for alignment characteristics ..."
echo -e "\nPlotting boxplots/histograms for alignment characteristics using alignmentSummary.R ...\n" >> R.log
if [ ! $location == "1" ]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f alignmentSummary.R
	R --slave -f alignmentSummary.R >> R.log 2>&1
else
	R --slave -f alignmentSummary.R >> R.log 2>&1
fi

#Plot boxplots/histograms for selected tree properties
cp gene_properties_update.txt tree_stats_table.csv
cat tree_stats_table.csv | awk '{ print $34 "," $35 "," $36 "," $37 "," $38 "," $39 "," $40 "," $41 }' > tmp && mv tmp tree_stats_table.csv

sed -i.bak 's/Bootstrap/Average_bootstrap/' tree_stats_table.csv
sed -i.bak 's/Branch_length/Average_branch_length/' tree_stats_table.csv
sed -i.bak 's/P_distance/Avg_p_dist/' tree_stats_table.csv
sed -i.bak 's/Satur_slope/Slope/' tree_stats_table.csv
sed -i.bak 's/Satur_R_sq/R_squared/' tree_stats_table.csv

echo -e "\nPlotting boxplots/histograms for tree properties..."
echo -e "\nPlotting boxplots/histograms for tree properties using treepropsPlot.r\n" >> R.log
if [ ! $location == "1" ]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f treepropsPlot.r
	R --slave -f treepropsPlot.r >> R.log 2>&1
else
	R --slave -f treepropsPlot.r >> R.log 2>&1
fi

#Copy all resulting PNGs to home
cp *histogram.png $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/
#Copy R.log to home
cp R.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/
#Prepare list of genes of updated selection
cat gene_properties_update.txt | sed 1d | cut -f1 | sort | sed 's/Corrected//g' | sed "s/_modif${MISSINGPERCENT}.fas//g" > selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt
mkdir -p $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes
cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes
echo -e "\nList of updated selected genes saved to $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt..."
echo -e "Change 'update' to 'update=yes' in 'settings.cfg' and continue with running scripts 7 and 8..."

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir09
fi

echo -e "\nScript HybPhyloMaker9 finished...\n"
