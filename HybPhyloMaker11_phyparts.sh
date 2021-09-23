#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=2:mem=16gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker11_PhyParts_Astral
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker11_PhyParts_Astral
#$ -o HybPhyloMaker11_PhyParts_Astral.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                 Script 11 - PhyParts for Astral species tree                 *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Runs PhyParts and PhyParts PieCharts for ASTRAL tree
#Takes all selected gene trees and only select trees with outgroup

#requires run on MetaCentrum or installation with 'install_software.sh'
#(i.e., 'phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar' must be in HybSeqSource (or in /usr/local/bin/)

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker11_PhyParts is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add phyparts-0.0.1
	module add newick-utils-13042016
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker11_PhyParts is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir11
	cd workdir11
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #adds NewickUtilities
	module load ???phyparts#???
	echo "PhyParts not yet installed on Hydra, exiting..."
	exit 3
else
	echo -e "\nHybPhyloMaker11_PhyParts is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir11
	cd workdir11
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
		echo -e "...with corrected reading frame"
	else
		mafftpath=$type/60mafft
		alnpath=$type/70concatenated_exon_alignments
		alnpathselected=$type/71selected
		treepath=$type/72trees
		echo -e ""
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/$selection/61mafft_corrected
		alnpath=$type/$selection/80concatenated_exon_alignments_corrected
		alnpathselected=$type/$selection/81selected_corrected
		treepath=$type/$selection/82trees_corrected
		echo -e "...with corrected reading frame...and for selection: $selection"
	else
		mafftpath=$type/$selection/60mafft
		alnpath=$type/$selection/70concatenated_exon_alignments
		alnpathselected=$type/$selection/71selected
		treepath=$type/$selection/72trees
		echo -e "...and for selection: $selection"
	fi
fi

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection"
fi

if [[ ! $collapse -eq "0" ]]; then
	echo -e "...and with trees with branches below ${collapse} BS collapsed"
else
	if [[ $requisite =~ "no" ]]; then
		echo -e "\n"
	fi
fi

if [[ $requisite =~ "yes" ]]; then
	echo -e "...and only with trees with requisite taxa present\n"
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir11)" ]; then
		echo -e "Directory 'workdir11' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir11 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM11
echo -e "HybPhyloMaker11: PhyParts for Astral species tree" > ${logname}.log
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "run on MetaCentrum: $PBS_O_HOST" >> ${logname}.log
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
ppcolors=${ppcolors// /-} #change ' ' in ${ppcolors} to '-'
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree OUTGROUP collapse requisite phypartsbs ppcolors nrpptrees; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
ppcolors=${ppcolors//-/ } #change '-' in ${ppcolors} back to ' '
if [[ $requisite =~ "yes" ]]; then
	echo -e "\nList of requisite samples" >> ${logname}.log
	echo $requisitetaxa | tr '|' '\n' >> ${logname}.log
fi
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Settings for collapsed and requisite selection
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		modif=with_requisite/collapsed${collapse}/
		treefile=trees_with_requisite_collapsed${collapse}.newick
	else
		modif=with_requisite/
		if [ -z "$OUTGROUP" ]; then
			treefile=trees_with_requisite.newick
		else
			treefile=trees_rooted_with_requisite.newick
		fi
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		modif=collapsed${collapse}/
		treefile=trees_collapsed${collapse}.newick
	else
		modif=""
		if [ -z "$OUTGROUP" ]; then
			treefile=trees.newick
		else
			treefile=trees_rooted.newick
		fi
	fi
fi

#Check if $OUTGROUP is set
if [ -z $OUTGROUP ]; then
	echo -e "Outgroup is not set. Exiting..."
	rm -d ../workdir11 2>/dev/null
	exit 3
fi

#Define Astral tree name
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite_collapsed${collapse}.tre
	else
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.tre
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}_collapsed${collapse}.tre
	else
		astraltree=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
	fi
fi

#Remove '.tre' from species tree name
sptree=$(cut -d'.' -f1 <<< $astraltree)

#Copy species tree
echo -e "\nCopying trees..."
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/${astraltree} .
else
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/${astraltree} .
fi

#Copy gene trees (file with all gene trees)
if [[ $update =~ "yes" ]]; then
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	fi
else
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	fi
fi

#Modify Astral tree (replace ' ' back to '-' and '_')
sed -i 's/ \([^ ]*\) / \1_/g' ${sptree}.tre #replace every second occurrence of ' ' by '_'
sed -i 's/ /-/g' ${sptree}.tre #replace all spaces by '-'

#Modify labels in gene tree
sed -i.bak 's/XX/-/g' $treefile
sed -i.bak2 's/YY/_/g' $treefile

#Get only rooted gene trees (i.e., trees containing outgroup)
nrgenetreesorig=$(wc -l < $treefile)
grep "$OUTGROUP" ${treefile} > tmp
mv tmp ${treefile}

#Check if there are any gene trees left
if [ $(wc -l < ${treefile}) -eq 0 ]; then
	echo -e "\nThere are no gene trees rooted with ${OUTGROUP}. Exiting..."
	exit 3
fi

#Reroot genestrees with $OUTGROUP
nw_reroot -s ${treefile} $OUTGROUP > tmp && mv tmp ${treefile}

#Put single tree per file to directory 'trees'
mkdir trees
split -a 4 -d -l 1 ${treefile} trees/tree_

nrgenetreesrooted=$(ls trees/tree_* | wc -l)

#Subselect only first 'nrpptrees' if necessary (to decrease running time)
if [ -z "$nrpptrees" ]; then #test whether $nrpptrees is empty
	nrpptrees=$nrgenetreesrooted #if empty set to number of rooted gene trees
fi

if [ $nrpptrees -gt $nrgenetreesrooted ]; then #test whether $nrpptrees is higher than nr. of rooted trees
	nrpptrees=$nrgenetreesrooted #if higher set to number of rooted gene trees
fi

if [ $nrpptrees -lt $nrgenetreesrooted ]; then
	echo -e "\nSubselecting $nrpptrees trees..."
	mkdir trees$nrpptrees
	cd trees
	ls tree_* | head -n $nrpptrees | xargs -I{} cp "{}" ../trees$nrpptrees/
	cd ..
	rm -r trees
	mv trees$nrpptrees trees
fi

#Calculate number of gene trees
nrgenetrees=$(ls trees/tree_* | wc -l)

#Set to 0.5 if $phypartsbs is empty
if [ -z "$phypartsbs" ]; then #test whether $phypartsbs is empty
	echo -e "\n'phypartsbs' is not set to any value, using 0.5..."
	phypartsbs=0.5
fi

#Statisctics
echo -e "\nSpecies tree: ${sptree}.tre" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Gene tree file: ${treefile}" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Nr. gene trees: ${nrgenetreesorig}" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Nr. gene trees rooted with '$OUTGROUP': ${nrgenetreesrooted}" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Using $nrpptrees gene trees" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt

#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
fi

#Run phyparts
echo -e "\nRunning PhyParts with support cutoff ${phypartsbs}..."
# -a what kind of analysis (0 - concon, 1 - fullconcon, 2 - duplications)
# -d directory of trees
# -m  mapping tree (species tree)
# -o prepend output files with this
# -s support cutoff (only keep things with greater support than the one specified)
# -v include verbose output
if [[ $location =~ "0" ]]; then
	#download java PhyParts
	cp /usr/local/bin/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar .
	cp $source/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar .
	echo "java -jar phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v"
	java -jar phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v > phyparts.log 2>&1
else
	echo "phyparts -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v"
	phyparts -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v > phyparts.log 2>&1
	#java -jar /software/phyparts/0.0.1/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v > phyparts.log 2>&1
fi

#Run phypartspiecharts
echo -e "\nRunning PhyParts PieCharts..."
echo -e "Using colours: ${ppcolors}"
# --svg_name name of the resulting graphics
# species tree
# prefix of phypart results (i.e., the same as '-o' option in phyparts)
# number of genetrees
if [[ $location =~ "0" ]]; then
	#download python script
	wget https://raw.githubusercontent.com/mossmatters/phyloscripts/master/phypartspiecharts/phypartspiecharts.py 2>/dev/null
	export QT_QPA_PLATFORM='offscreen'
	if [ -z "$ppcolors" ]; then
		echo "python3 phypartspiecharts.py --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees}"
		python3 phypartspiecharts.py --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} 2>/dev/null
	else
		echo "python3 phypartspiecharts.py --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}"
		python3 phypartspiecharts.py --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors} 2>/dev/null
		#python /software/phyparts/0.0.1/target/phypartspiecharts.py --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}
	fi
else
	if [ -z "$ppcolors" ]; then
		echo "phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees}"
		phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees}
	else
		echo "phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}"
		phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}
		#python /software/phyparts/0.0.1/target/phypartspiecharts.py --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}
	fi
fi

#Convert SVG to PDF
if [[ $PBS_O_HOST == *".cz" ]]; then
	module add python36-modules-gcc #adds also cairosvg
fi
cairosvg phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg -o phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.pdf

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}
	cp phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}
	cp phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}
	cp trees_res.*node* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
	cp trees_res.*hist* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
	cp phyparts.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
else
	cp phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}
	cp phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.pdf $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}
	cp phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}
	cp trees_res.*node* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
	cp trees_res.*hist* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
	cp phyparts.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}Astral/phyparts_${phypartsbs}
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}Astral/phyparts_${phypartsbs}
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir11
fi

echo -e "\nScript HybPhyloMaker11 finished...\n"
