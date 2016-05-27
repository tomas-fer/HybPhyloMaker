#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -N HybPipe9_update_trees
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe9_update_trees
#$ -o HybPipe9_update_trees.log

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                           Script 09 - Update trees                           *
# *                                   v.1.0.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#UPDATE on tree selection
# 

if [[ $PBS_O_HOST == *".cz" ]]; then
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
	echo "Hydra..."
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
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir09
	cd workdir09
fi
#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=80concatenated_exon_alignments_corrected
	alnpathselected=81selected_corrected
	treepath=82trees_corrected
else
	alnpath=70concatenated_exon_alignments
	alnpathselected=71selected
	treepath=72trees
fi

#Copy scripts
cp $source/plotting_correlations.R .
cp $source/alignmentSummary.R .
cp $source/treepropsPlot.r .

#Copy updated gene list with properties (with unwanted genes deleted)
cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/gene_properties_update.txt .
#Plot gene properties correlations
cp gene_properties_update.txt combined.txt
echo -e "Plotting gene properties correlations for updated selection...\n"
if [ ! $LOGNAME == "" ]; then
	#Run R script for correlation visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#echo "Running xvfb-run R..."
	#xvfb-run R --slave -f plotting_correlations.R
	R --slave -f plotting_correlations.R
else
	R --slave -f plotting_correlations.R
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
if [ ! $LOGNAME == "" ]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f alignmentSummary.R
	R --slave -f alignmentSummary.R
else
	R --slave -f alignmentSummary.R
fi

#Plot boxplots/histograms for selected tree properties
cp gene_properties_update.txt tree_stats_table.csv
cat tree_stats_table.csv | awk '{ print $34 "," $35 "," $36 "," $37 "," $38 "," $39 "," $40 "," $41 }' > tmp && mv tmp tree_stats_table.csv

sed -i.bak 's/Bootstrap/Average_bootstrap/' tree_stats_table.csv
sed -i.bak 's/Branch_length/Average_branch_length/' tree_stats_table.csv
sed -i.bak 's/P_distance/Avg_p_dist/' tree_stats_table.csv
sed -i.bak 's/Satur_slope/Slope/' tree_stats_table.csv
sed -i.bak 's/Satur_R_sq/R_squared/' tree_stats_table.csv

echo -e "Plotting boxplots/histograms for tree properties...\n"
if [ ! $LOGNAME == "" ]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f treepropsPlot.r
	R --slave -f treepropsPlot.r
else
	R --slave -f treepropsPlot.r
fi

#Copy all resulting PNGs to home
cp *histogram.png $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/

#Prepare list of genes of updated selection
cat gene_properties_update.txt | sed 1d | cut -f1 | sort | sed 's/Corrected/CorrectedAssembly_/g' | sed 's/_modif70.fas//g' > selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt
mkdir $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes
cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir09
fi
