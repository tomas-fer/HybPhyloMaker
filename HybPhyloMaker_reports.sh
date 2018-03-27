#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=0:10:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker_reports
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 2
#$ -q sThM.q
#$ -l mres=1G,h_data=1G,h_vmem=1G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker_reports
#$ -o HybPhyloMaker_reports.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *          Additional script - combine summary tables into single XLS          *
# *                                   v.1.6.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2018 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#This script combines summary tables into a single XLS file
# - filtering summary (reads_summary.txt)
# - mapping summary (mapping_summary.txt)
# - missing data overview for all genes (MissingDataOverview.txt)
# - missing data overview for selected genes (MissingDataOverview_${MISSINGPERCENT}.txt)
# - alignment/tree characteristics for all genes (summaryALL.txt)
# - alignment/tree characteristics for selected genes (gene_properties.txt)
#Run this after building gene trees with 6 (and calculating properties with 6a2 in case of RAxML)

#Requires 'openxlsx' R package and HybPhyloMaker_reports.R (in 'HybSeqSource' folder)

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker-summary is running on MetaCentrum..."
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
	module add R-3.4.3-gcc
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker-summary is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdirsummary
	cd workdirsummary
	#Add necessary modules
	module load tools/R/3.3.1
else
	echo -e "\nHybPhyloMaker-summary is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdirsummary
	cd workdirsummary
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
	echo -e "...and with updated gene selection"
else
	if [[ $requisite =~ "no" ]]; then
		echo -e ""
	fi
fi

if [[ $requisite =~ "yes" ]]; then
	echo -e "...and only with trees with requisite taxa present\n"
else
	echo -e "\n"
fi

#Copy R script
cp $source/HybPhyloMaker_reports.R .

#Copy all tables
echo "Collecting tables..."
cp ${path}/20filtered/reads_summary.txt .
cp ${path}/${type}/21mapped_${mappingmethod}/mapping_summary.txt .
cp ${path}/${alnpathselected}${MISSINGPERCENT}/MissingDataOverview.txt .
cp ${path}/${alnpathselected}${MISSINGPERCENT}/MissingDataOverview_${MISSINGPERCENT}.txt .
cp ${path}/${alnpathselected}${MISSINGPERCENT}/summaryALL.txt .
cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/gene_properties.txt .

#Change to csv
echo "Modifying tables..."
sed "s/\t/,/g" reads_summary.txt > 1-Reads_summary.csv #Change TABs to commas
sed "s/\t/,/g" mapping_summary.txt > 2-Mapping_summary.csv #Change TABs to commas

#Transpose MissingDataOverview
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
}' MissingDataOverview.txt > tmp && mv tmp MissingDataOverview.txt

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
}' MissingDataOverview_${MISSINGPERCENT}.txt > tmp && mv tmp MissingDataOverview_${MISSINGPERCENT}.txt

#Change to csv and modify
sed "s/ /,/g" MissingDataOverview.txt > 3-Missing_data.csv #Change spaces to commas
sed -i.bak '1s/-/x/g' 3-Missing_data.csv #Change '-' to 'x' (only on the first line, i.e., in gene names)
sed -i.bak2 's/Assembly_//g' 3-Missing_data.csv
sed -i.bak3 's/_/ /g' 3-Missing_data.csv

sed "s/ /,/g" MissingDataOverview_${MISSINGPERCENT}.txt > 4-Missing_data_${MISSINGPERCENT}.csv #Change spaces to commas
sed -i.bak '1s/-/x/g' 4-Missing_data_${MISSINGPERCENT}.csv #Change '-' to 'x' (only on the first line, i.e., in gene names)
sed -i.bak2 's/Assembly_//g' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak3 's/_/ /g' 4-Missing_data_${MISSINGPERCENT}.csv

sed "s/\t/,/g" summaryALL.txt > 5-All_genes.csv #Change TABs to commas
sed -i.bak 's/Assembly_//g' 5-All_genes.csv
sed -i.bak2 's/\.fasta//g' 5-All_genes.csv
sed -i.bak3 's/_/ /g' 5-All_genes.csv
sed -i.bak4 's/inf,/NA,/g' 5-All_genes.csv #Change 'inf' (introduced by MstatX) to 'NA'

sed "s/\t/,/g" gene_properties.txt > 6-Selected_genes.csv #Change TABs to commas
sed -i.bak 's/\.fas//g' 6-Selected_genes.csv
sed -i.bak2 's/_/ /g' 6-Selected_genes.csv

#Modify columns
#3-Missing_data.csv
sed -i.bak4 's/average missing data/,,averagemissingdata/' 3-Missing_data.csv
sed -i.bak5 's/nr assemblies with completely missing data/,,nrassemblieswithcompletelymissingdata/' 3-Missing_data.csv
sed -i.bak6 's/ /,/g' 3-Missing_data.csv #Change spaces to commas
sed -i.bak7 's/-/,/g' 3-Missing_data.csv #Change '-' (between genus-species) to commas (in order to make separate column)
head=$(head -1 3-Missing_data.csv) #First line to $head
echo "Sample no.,Genus,$head" > head.txt #Modify $head (add two columns)
sed -i.bak 's/species/Species/' head.txt
sed -i.bak8 '1d' 3-Missing_data.csv #Delete first line
cat 3-Missing_data.csv | cut -d"," -f1,2 > b.txt #First and second columns (genus, species)
cat 3-Missing_data.csv | cut -d"," -f3 > a.txt #Third column (code)
cat 3-Missing_data.csv | cut -d"," -f4- > c.txt #Fourth and remaining columns (data)
paste a.txt b.txt c.txt > comb.txt #Combine columns back
cat head.txt comb.txt > 3-Missing_data.csv #Combine with head
sed -i.bak9 "s/\t/,/g" 3-Missing_data.csv #Change TABs to commas
sed -i.bak10 's/averagemissingdata/Average missing data/' 3-Missing_data.csv
sed -i.bak11 's/nrassemblieswithcompletelymissingdata/Nr. genes with completely missing data/' 3-Missing_data.csv
rm a.txt b.txt c.txt head.txt comb.txt *.bak*

#4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak4 's/averageMissing/,,averageMissing/' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak5 's/percPresentSpecies/percPresentSpecies,,/' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak6 's/average missing/averageMissing/' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak7 's/total genes/totalgenes/' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak8 's/ /,/g' 4-Missing_data_${MISSINGPERCENT}.csv #Change spaces to commas
tail -1 4-Missing_data_${MISSINGPERCENT}.csv > tail.txt #Last line to tail.txt
sed -i.bak9 '$ d' 4-Missing_data_${MISSINGPERCENT}.csv #Delete last line
sed -i.bak10 's/-/,/g' 4-Missing_data_${MISSINGPERCENT}.csv #Change '-' (between genus-species) to commas (in order to make separate column)
head=$(head -1 4-Missing_data_${MISSINGPERCENT}.csv) #First line to $head
echo "Sample no.,Genus,$head" > head.txt #Modify $head (add two columns)
sed -i.bak 's/species/Species/' head.txt
sed -i.bak11 '1d' 4-Missing_data_${MISSINGPERCENT}.csv #Delete first line
cat 4-Missing_data_${MISSINGPERCENT}.csv | cut -d"," -f1,2 > b.txt #First and second columns (genus, species)
cat 4-Missing_data_${MISSINGPERCENT}.csv | cut -d"," -f3 > a.txt #Third column (code)
cat 4-Missing_data_${MISSINGPERCENT}.csv | cut -d"," -f4- > c.txt #Fourth and remaining columns (data)
paste a.txt b.txt c.txt > comb.txt #Combine columns back
cat head.txt comb.txt tail.txt> 4-Missing_data_${MISSINGPERCENT}.csv #Combine with head and tail
sed -i.bak12 "s/\t/,/g" 4-Missing_data_${MISSINGPERCENT}.csv #Change TABs to commas
sed -i.bak13 's/averageMissing/Average missing data/g' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak14 's/percPresentSpecies/Species presence (%)/' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak15 's/totalgenes/Total nr. of genes/' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i.bak16 's/N\/A/NA/g' 4-Missing_data_${MISSINGPERCENT}.csv #Change 'N/A' to 'NA' (to be recognized by R as missing data)
rm a.txt b.txt c.txt head.txt comb.txt tail.txt *.bak*

#Get number of samples, genes and selected genes
nrsamples=$(cat $path/10rawreads/SamplesFileNames.txt | wc -l)
if [[ $cp =~ "yes" ]]; then
	nrgenes=$(ls ${path}/cp/60mafft/*.fasta | wc -l)
else
	nrgenes=$(ls ${path}/${alnpath}/*.fasta | wc -l)
fi
nrselected=$(cat ${path}/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | wc -l)

#Make XLS file
echo "Writing XLSX..."
# This requires pyexcel-cli and pyexcel-xlsx (to be installed with pip install)
# This produces unformated XLSX
# `python <<END
# from pyexcel.cookbook import merge_all_to_a_book;
# import glob;
# merge_all_to_a_book(glob.glob("*.csv"), "summary.xlsx");
# END`
R --slave -f HybPhyloMaker_reports.R $nrsamples $nrgenes $nrselected $MISSINGPERCENT

mv summary.xlsx ${data}_summary.xlsx

#Copy summary to home
cp ${data}_summary.xlsx ${path}

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdirsummary
fi

echo -e "\nScript HybPhyloMaker-summary finished...\n"
