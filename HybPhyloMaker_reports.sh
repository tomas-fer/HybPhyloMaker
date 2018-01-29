#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=4:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker10_requisite_collapse
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 2
#$ -q sThM.q
#$ -l mres=8G,h_data=8G,h_vmem=8G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker10_requisite_collapse
#$ -o HybPhyloMaker10_requisite_collapse.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *          Additional script - combine summary tabels into single XLS          *
# *                                   v.1.5.1                                    *
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
sed "s/\t/,/g" reads_summary.txt > 1-Reads_summary.csv
sed "s/\t/,/g" mapping_summary.txt > 2-Mapping_summary.csv

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
sed "s/ /,/g" MissingDataOverview.txt > 3-Missing_data.csv
sed -i 's/Assembly_//g' 3-Missing_data.csv
sed -i 's/_/ /g' 3-Missing_data.csv
sed "s/ /,/g" MissingDataOverview_${MISSINGPERCENT}.txt > 4-Missing_data_${MISSINGPERCENT}.csv
sed -i 's/Assembly_//g' 4-Missing_data_${MISSINGPERCENT}.csv
sed -i 's/_/ /g' 4-Missing_data_${MISSINGPERCENT}.csv

sed "s/\t/,/g" summaryALL.txt > 5-All_genes.csv
sed -i 's/Assembly_//g' 5-All_genes.csv
sed -i 's/\.fasta//g' 5-All_genes.csv
sed -i 's/_/ /g' 5-All_genes.csv
sed "s/\t/,/g" gene_properties.txt > 6-Selected_genes.csv
sed -i 's/\.fas//g' 6-Selected_genes.csv
sed -i 's/_/ /g' 6-Selected_genes.csv

#Make XLS file
echo "Writing XLSX..."
`python <<END
from pyexcel.cookbook import merge_all_to_a_book;
import glob;
merge_all_to_a_book(glob.glob("*.csv"), "Summary.xlsx");
END`

mv Summary.xlsx ${data}_summary.xlsx

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
