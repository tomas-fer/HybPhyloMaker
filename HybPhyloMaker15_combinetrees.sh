#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker15_combinetrees
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker15_combinetrees
#$ -o HybPhyloMaker15_combinetrees.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                Script 15 - Combine trees with support values                 *
# *                                   v.1.8.0b                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Combines bootstrap support values from two trees into a single tree...
#Now works for Astral, Astral4, Astrid, MRL, FastTree and ExaML

if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker15 is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	module add py-p4phylogenetics/20240606
	module add newick-utils-13042016
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker15 is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir15
	cd workdir15
else
	echo -e "\nHybPhyloMaker15 is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir15
	cd workdir15
fi

if [[ $cp =~ "yes" ]]; then
	echo -e "This script does not work with cp data yet. Exiting..."
	rm -d ../workdir15 2>/dev/null
	exit 3
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir15)" ]; then
		echo -e "Directory 'workdir15' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir15 2>/dev/null
		exit 3
	fi
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "\nWorking with cpDNA"
	type="cp"
else
	echo -en "\nWorking with exons"
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

if [[ $update =~ "yes" ]]; then
	echo -e "...and with updated gene selection"
else
	echo -e ""
fi

if [[ ! $collapse -eq "0" ]]; then
	echo -e "...and with trees with branches below ${collapse} BS collapsed"
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

#Settings for collapsed and requisite selection
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		modif=with_requisite/collapsed${collapse}/
		modif1=with_requisite/collapsed${collapse}
		treefile=trees_with_requisite_collapsed${collapse}.newick
	else
		modif=with_requisite/
		modif1=with_requisite
		if [ -z "$OUTGROUP" ]; then
			treefile=trees_with_requisite.newick
		else
			treefile=trees_rooted_with_requisite.newick
		fi
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		modif=collapsed${collapse}/
		modif1=collapsed${collapse}
		treefile=trees_collapsed${collapse}.newick
	else
		modif=""
		modif1=""
		if [ -z "$OUTGROUP" ]; then
			treefile=trees.newick
		else
			treefile=trees_rooted.newick
		fi
	fi
fi

#Write log
logname=HPM15
echo -e "HybPhyloMaker15: Combine species trees with support values" > ${logname}.log
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "Job run on MetaCentrum: $PBS_JOBID" >> ${logname}.log
	echo -e "From: $PBS_O_HOST" >> ${logname}.log
	echo -e "Host: $HOSTNAME" >> ${logname}.log
	echo -e "$PBS_NUM_NODES node(s) with $PBS_NCPUS core(s)" >> ${logname}.log
	memM=$(bc <<< "scale=2; $(echo $PBS_RESC_MEM) / 1024 / 1024 ")
	memG=$(bc <<< "scale=2; $(echo $PBS_RESC_MEM) / 1024 / 1024 / 1024 ")
	if (( $(echo $memG 1 | awk '{if ($1 < $2) print 1;}') )); then
		echo -e "Memory: $memM Mb" >> ${logname}.log
	else
		echo -e "Memory: $memG Gb" >> ${logname}.log
	fi
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
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree OUTGROUP requisite tree1 tree2 prec; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [[ $requisite =~ "yes" ]]; then
	echo -e "\nList of requisite samples" >> ${logname}.log
	echo $requisitetaxa | tr '|' '\n' >> ${logname}.log
fi
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Copy the python3 script
cp /storage/$server/home/$LOGNAME/HybSeqSource/combineboot3.py .

#settings based on tree selection
if [[ $tree1 = "Astral" ]]; then
	tpath1=Astral
	tname1=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
elif [[ $tree1 = "Astral4" ]]; then
	tpath1=Astral4
	tname1=Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
elif [[ $tree1 = "FastTree" ]]; then
	tpath1=concatenated
	tname1=concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
elif [[ $tree1 = "ExaML" ]]; then
	tpath1=concatenatedExaML
	tname1=ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre
elif [[ $tree1 = "MRL" ]]; then
	tpath1=MRL
	tname1=MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
elif [[ $tree1 = "Astrid" ]]; then
	tpath1=Astrid
	tname1=Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
else
	echo -e "Variable 'tree1' is not set to one of allowed values (Astral, Astral4, Astrid, MRL, FastTree, ExaML). Exiting...\n"
	rm -d ../workdir15/ 2>/dev/null
	exit 3
fi

if [[ $tree2 = "Astral" ]]; then
	tpath2=Astral
	tname2=Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
elif [[ $tree2 = "Astral4" ]]; then
	tpath2=Astral4
	tname2=Astral4_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
elif [[ $tree2 = "FastTree" ]]; then
	tpath2=concatenated
	tname2=concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
elif [[ $tree2 = "ExaML" ]]; then
	tpath2=concatenatedExaML
	tname2=ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre
elif [[ $tree2 = "MRL" ]]; then
	tpath2=MRL
	tname2=MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
elif [[ $tree2 = "Astrid" ]]; then
	tpath2=Astrid
	tname2=Astrid_${MISSINGPERCENT}_${SPECIESPRESENCE}_withbootstrap.tre
else
	echo -e "Variable 'tree2' is not set to one of allowed values (Astral, Astral4, Astrid, MRL, FastTree, ExaML). Exiting...\n"
	rm -d ../workdir15/ 2>/dev/null
	exit 3
fi

#Copy the two trees
echo -e "Copying trees..."
echo -e " tree 1:\t$tname1"
echo -e " tree 2:\t$tname2"
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpath1}/${tname1} .
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${tpath2}/${tname2} .
else
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpath1}/${tname1} .
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${tpath2}/${tname2} .
fi

#Modify names in trees (in order to have 'GenusXXspeciesXXcode')
echo -e "\nModifying names..."
mv ${tname1} tree1.tre
mv ${tname2} tree2.tre

if [[ $tree1 =~ "ExaML" ]]; then
	sed -i 's/-/XX/g' tree1.tre #change the names in tree1
	sed -i 's/_/XX/g' tree1.tre #change the names in tree1
else
	sed -i 's/ /XX/g' tree1.tre #change the names in tree1
fi

if [[ $tree2 =~ "ExaML" ]]; then
	sed -i 's/-/XX/g' tree2.tre #change the names in tree2
	sed -i 's/_/XX/g' tree2.tre #change the names in tree2
else
	sed -i 's/ /XX/g' tree2.tre #change the names in tree2
fi

#Reroot trees
#modify ${OUTGROUP} first ('-' and '_' to 'XX')
OUTGROUP=${OUTGROUP/-/XX}
OUTGROUP=${OUTGROUP/_/XX}
nw_reroot -s tree1.tre $OUTGROUP > tmp && mv tmp tree1.tre
nw_reroot -s tree2.tre $OUTGROUP > tmp && mv tmp tree2.tre

#Round support values in Astral or Astral4 to '$prec' decimals only
echo -e "\nRounding support values to $prec decimals..."
if [[ $tree1 =~ "Astral" || $tree1 =~ "Astral4" || $tree1 =~ "FastTree" ]]; then
	nw_labels -L tree1.tre > bs #export BS values
	#this does not correct rounding
	#nw_labels -L tree1.tre | numfmt --field=1- --format=%.${prec}f | grep -o '.*[1-9]' > bsr #export BS values, round to ${prec} decimals and trim off trailing zeroes (grep command)
	#AWK is better ('*' before 'f' takes the value after ',')
	nw_labels -L tree1.tre | awk -v a=${prec} '{printf "%.*f\n", a, $1}' | grep -o '.*[1-9]' > bsr #export BS values, round to ${prec} decimals and trim off trailing zeroes (grep command)
	paste bs bsr > replace #combine the two files
	#replace BS in tree1 by rounded BS
	cat replace | while read -r a b; do
		sed -i "s/$a/$b/" tree1.tre
	done
fi

if [[ $tree2 =~ "Astral" || $tree2 =~ "Astral4" || $tree2 =~ "FastTree" ]]; then
	nw_labels -L tree2.tre > bs2 #export BS values
	#this does not correct rounding
	#nw_labels -L tree1.tre | numfmt --field=1- --format=%.${prec}f | grep -o '.*[1-9]' > bsr #export BS values, round to ${prec} decimals and trim off trailing zeroes (grep command)
	#AWK is better ('*' before 'f' takes the value after ',')
	nw_labels -L tree2.tre | awk -v a=${prec} '{printf "%.*f\n", a, $1}' | grep -o '.*[1-9]' > bsr2 #export BS values, round to ${prec} decimals and trim off trailing zeroes (grep command)
	paste bs2 bsr2 > replace2 #combine the two files
	#replace BS in tree2 by rounded BS
	cat replace2 | while read -r a b; do
		sed -i "s/$a/$b/" tree2.tre
	done
fi

#Copy & rename & modify alignment in phylip format
echo -e "\nCopying alignment..."
if [[ $update =~ "yes" ]]; then
	if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
	else
		echo -e "No concatenated alignment file in PHYLIP format found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenated'. Exiting..."
		rm -d ../workdir15 2>/dev/null
		exit 3
	fi
else
	if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
	else
		echo -e "No concatenated alignment file in PHYLIP format found in '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenated'. Exiting..."
		rm -d ../workdir15 2>/dev/null
		exit 3
	fi
fi
mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip concatenated.phylip #the name has to be concatenated.phylip
sed -i 's/_contigs.fas//' concatenated.phylip #change the names in alignment (replace '_contigs.fas' by nothing)
sed -i 's/-/XX/' concatenated.phylip #change the names in alignment (replace first occurrence of '-' at each line by 'XX')
sed -i 's/_/XX/' concatenated.phylip #change the names in alignment

#Run the python script
echo -e "\nRunning combination of support values..."
python3 ./combineboot3.py tree1.tre tree2.tre
#it creates 'combinedSupportsTree.nex' (NEXUS format)

#Transform the tree to NEWICK and modify
grep "&" combinedSupportsTree.nex | cut -d" " -f7- > combinedSupportsTree.tre
sed -i 's/ //g' combinedSupportsTree.tre #remove all gaps
sed -i 's/XX/ /g' combinedSupportsTree.tre #XX to gaps
sed -i 's/ \([^ ]*\) / \1_/g' combinedSupportsTree.tre #every second gap to '_'
sed -i 's/ /-/g' combinedSupportsTree.tre #remaining gaps to '-'

#Replace '1.000' to '1'
sed -i 's/1.000/1/g' combinedSupportsTree.tre

#Replace full supports by '*', e.g. '1/100' to '*/*'
sed -e "s/\/100'/\/*'/g" -e "s/'100\//'*\//g" -e "s/'1\//'*\//g" -e "s/\/1'/\/*'/g" combinedSupportsTree.tre > ${tree1}_and_${tree2}_${MISSINGPERCENT}_${SPECIESPRESENCE}_stars.tre
#Remove full supports, i.e. '*/*'
sed "s/\*\/\*//g" ${tree1}_and_${tree2}_${MISSINGPERCENT}_${SPECIESPRESENCE}_stars.tre > ${tree1}_and_${tree2}_${MISSINGPERCENT}_${SPECIESPRESENCE}_noStars.tre
#Rename the final combined tree
mv combinedSupportsTree.tre ${tree1}_and_${tree2}_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
#Copy the final tree to home
if [[ $update =~ "yes" ]]; then
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}combinedTrees
	cp ${tree1}_and_${tree2}_${MISSINGPERCENT}_${SPECIESPRESENCE}*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}combinedTrees
else
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}combinedTrees
	cp ${tree1}_and_${tree2}_${MISSINGPERCENT}_${SPECIESPRESENCE}*.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}combinedTrees
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}combinedTrees
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}combinedTrees
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir15
fi

echo -e "\nScript HybPhyloMaker 15 finished...\n"

