#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker8m2_PhyloNet_summary
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=2G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8m2_PhyloNet
#$ -o HybPhyloMaker8m2_PhyloNet.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                        Script 08m2 - PhyloNet summary                        *
# *                                   v.1.8.0                                    *
# *                                  Tomas Fer                                   *
# * Dept. of Botany, Charles University, Prague, Czech Republic, 2023            *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#This summarizes finished PhyloNet runs, creates summary table and plots resulting networks

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8m2 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add r/4.0.0-gcc
	module add julia
	module add newick-utils-13042016
	module add python36-modules-gcc #to add cairosvg
	module add jdk-8
	#module add debian9-compat
	#set nr processors
	#cpu=$TORQUE_RESC_PROC
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8m2 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08m2
	cd workdir08m2
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #adds NewickUtilities
	module load tools/R/3.4.1
	module load java/1.7
	#module load bioinformatics/newickutilities/0.0
	#julia
	#cairosvg
	#cpu=$NSLOTS
else
	echo -e "\nHybPhyloMaker8m2 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08m2
	cd workdir08m2
	#cpu=4
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -en "Working with exons"
	type="exons"
fi

#Settings for selection and (un)corrected reading frame
if [ -z "$selection" ]; then
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/61mafft_corrected
		alnpath=$type/80concatenated_exon_alignments_corrected
		alnpathselected=$type/81selected_corrected
		treepath=$type/82trees_corrected
		echo -en "...with corrected reading frame"
	else
		mafftpath=$type/60mafft
		alnpath=$type/70concatenated_exon_alignments
		alnpathselected=$type/71selected
		treepath=$type/72trees
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/$selection/61mafft_corrected
		alnpath=$type/$selection/80concatenated_exon_alignments_corrected
		alnpathselected=$type/$selection/81selected_corrected
		treepath=$type/$selection/82trees_corrected
		echo -en "...with corrected reading frame...and for selection: $selection"
	else
		mafftpath=$type/$selection/60mafft
		alnpath=$type/$selection/70concatenated_exon_alignments
		alnpathselected=$type/$selection/71selected
		treepath=$type/$selection/72trees
		echo -en "...and for selection: $selection"
	fi
fi

#Write log
logname=HPM8m2
echo -e "HybPhyloMaker8m2: PhyloNet network summary" > ${logname}.log
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "run on MetaCentrum: $HOSTNAME (submitted from $PBS_O_HOST)" >> ${logname}.log
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
for set in data MISSINGPERCENT SPECIESPRESENCE selection cp corrected update tree hstart hmax; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

data1=$(echo $data | sed 's:.*/::') #everything after last slash, for naming purposes

#Check necessary file
echo -e "\n\nTesting if input data are available..."
if [[ $update =~ "yes" ]]; then
	for i in $(seq $hstart $hmax); do
		if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${i}/${data1}_${i}_reticulations.networks" ]; then
			echo -e "Results for ${i} reticulations OK"
		else
			echo -e "no networks for ${i} reticulations found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${i}'. Complete first PhyloNet run. Exiting..."
			rm -d ../workdir08m2 2>/dev/null
			exit 3
		fi
	done
	echo -e "\n"
else
	for i in $(seq $hstart $hmax); do
		if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${i}/${data1}_${i}_reticulations.networks" ]; then
			echo -e "Results for ${i} reticulations OK"
		else
			echo -e "no networks for ${i} reticulations found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${i}'. Complete first PhyloNet run. Exiting..."
			rm -d ../workdir08m2 2>/dev/null
			exit 3
		fi
	done
	echo -e "\n"
fi

#Copy script
cp $source/plotNetworks.jl .

#Summarize PhyloNet runs
for pn in $(seq $hstart $hmax); do
	#copy data from home
	if [[ $update =~ "yes" ]]; then
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}/${data1}_${pn}_reticulations.networks .
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}/${pn}_table.txt .
	else
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}/${data1}_${pn}_reticulations.networks .
		cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}/${pn}_table.txt .
	fi
	#put every network to a separate file
	grep "^(" ${data1}_${pn}_reticulations.networks | split -l 1 - -d -a 1 ${pn}reti_network
done

#Make a summary table
#combine result tables for all hybridizations
cat *_table.txt > ${data1}_PhyloNet_summary.txt
rm *_table.txt
#sum values in the first three columns and add it as the last column (ki values)
awk '{ sum = $1 + $2 + $3; print $0 "\t" sum }' ${data1}_PhyloNet_summary.txt > test && mv test ${data1}_PhyloNet_summary.txt
#calculate AIC=-2loglik + 2ki and add it as last column
#OFMT defines that values will be printed as decimals (not in scientific format), aic after a comma introduces a spaces that is then deleted with a sed command
awk -v OFMT='%f' '{ aic = -2 * $4 + 2 * $5 ; print $0 "\t",aic }' ${data1}_PhyloNet_summary.txt | sed 's/ //' > test && mv test ${data1}_PhyloNet_summary.txt
#calculate deltaAIC (subtract minimum AIC from AIC and add it as last column
minaic=$(awk '{if (min == "") min=$6 ; else if ($6 < min) min=$6}END{print min}' ${data1}_PhyloNet_summary.txt)
awk -v val=$minaic '{ delta = $6 - val ; print $0 "\t" delta}' ${data1}_PhyloNet_summary.txt > test && mv test ${data1}_PhyloNet_summary.txt
#sort according deltaAIC (i.e., column 7)
sort -n -k 7 ${data1}_PhyloNet_summary.txt > test && mv test ${data1}_PhyloNet_summary.txt
#add header
echo -e "NrHybridizations\tNrBranches\tNrGeneTrees\tlogLikelihood\tki\tAIC\tdeltaAIC" > header.txt
cat header.txt ${data1}_PhyloNet_summary.txt > test && mv test ${data1}_PhyloNet_summary.txt
rm header.txt

#Plot networks (in Julia)
julia plotNetworks.jl "$OUTGROUP"

#Create PDF from all SVG
echo -e "Transforming SVG to PDF...\n"
for i in $(ls *.svg | cut -d'.' -f1); do
	cairosvg ${i}.svg -o ${i}.pdf
done

#Merge PDFs and copy results to home
#download PDFbox
wget https://dlcdn.apache.org/pdfbox/2.0.30/pdfbox-app-2.0.30.jar
#loop over reticulations
for pn in $(seq $hstart $hmax); do
	java -jar pdfbox-app-2.0.30.jar PDFMerger ${pn}reti*_unrooted.pdf ${data1}_PhyloNet_${pn}_unrooted.pdf
	java -jar pdfbox-app-2.0.30.jar PDFMerger ${pn}reti*_rooted.pdf ${data1}_PhyloNet_${pn}_rooted.pdf
	if [[ $update =~ "yes" ]]; then
		cp ${pn}reti*.{pdf,svg} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}
		cp ${data1}_PhyloNet_${pn}_unrooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}
		cp ${data1}_PhyloNet_${pn}_rooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}
	else
		cp ${pn}reti*.{pdf,svg} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}
		cp ${data1}_PhyloNet_${pn}_unrooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}
		cp ${data1}_PhyloNet_${pn}_rooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}
	fi
done
#merge merged PDFs
java -jar pdfbox-app-2.0.30.jar PDFMerger ${data1}_PhyloNet_*_unrooted.pdf ${data1}_PhyloNet_unrooted.pdf
java -jar pdfbox-app-2.0.30.jar PDFMerger ${data1}_PhyloNet_*_rooted.pdf ${data1}_PhyloNet_rooted.pdf
#copy to home
if [[ $update =~ "yes" ]]; then
	cp ${data1}_PhyloNet_unrooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/
	cp ${data1}_PhyloNet_rooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/
	cp ${data1}_PhyloNet_summary.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/
else
	cp ${data1}_PhyloNet_unrooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/
	cp ${data1}_PhyloNet_rooted.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/
	cp ${data1}_PhyloNet_summary.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
	if [[ $update =~ "yes" ]]; then
		cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet
	else
		cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet
	fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08m2
	rm -r workdir08m2
fi

echo -e "\nHybPhyloMaker8m2 finished...\n"
