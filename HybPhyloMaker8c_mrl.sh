#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:0:0
#PBS -l select=1:ncpus=8:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8c_MRL
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8c_MRL
#$ -o HybPhyloMaker8c_MRL.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                        Script 08c - MRL species tree                         *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using MRL methods using RAxML from trees saved in single gene tree file (with *.newick suffix)
#RAxML options: standard/rapid bootstrap (raxmlboot=), number of bootstrap replicates (bsrep=)
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6a_RAxML_for_selected.sh or HybPhyloMaker6b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPhyloMaker7_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#Works also for trees after update, requisite taxa selection and collapsing (see HybPhyloMaker9_update_trees.sh and HybPhyloMaker10_requisite_collapse.sh)

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8c is running on MetaCentrum..."
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
	module add jdk-1.6.0
	module add raxml-8.2.4
	module add newick-utils-13042016
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8c is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08c
	cd workdir08c
	#Add necessary modules
	module load java/1.7
	module load bioinformatics/raxml/8.2.11
	module load bioinformatics/anaconda3/5.1 #NewickUtilities
else
	echo -e "\nHybPhyloMaker8c is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08c
	cd workdir08c
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

#Check necessary file
if [[ $requisite =~ "no" ]] && [[ $collapse -eq "0" ]];then
	echo -ne "Testing if input data are available..."
	if [[ $update =~ "yes" ]]; then
		if [ -z "$OUTGROUP" ]; then
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick" ]; then
				echo -e "OK\n"
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08c 2>/dev/null
				exit 3
			fi
		else
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
				echo -e "OK\n"
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08c 2>/dev/null
				exit 3
			fi
		fi
	else
		if [ -z "$OUTGROUP" ]; then
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick" ]; then
				echo -e "OK\n"
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08c 2>/dev/null
				exit 3
			fi
		else
			if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick" ]; then
				echo -e "OK\n"
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick' is missing. Exiting...\n"
				rm -d ../workdir08c 2>/dev/null
				exit 3
			fi
		fi
	fi
	
	#Test if folder for results exits
	if [[ $update =~ "yes" ]]; then
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/MRL" ]; then
			echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/MRL' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir08c 2>/dev/null
			exit 3
		fi
	else
		if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/MRL" ]; then
			echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/MRL' already exists. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir08c 2>/dev/null
			exit 3
		fi
	fi
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08c)" ]; then
		echo -e "Directory 'workdir08c' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08c 2>/dev/null
		exit 3
	fi
fi

#Add necessary programs and files
cp $source/mrp.jar .

#Copy genetree file
if [[ $update =~ "yes" ]]; then
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
		#mv trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	fi
else
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
		#mv trees${MISSINGPERCENT}_${SPECIESPRESENCE}_withoutBS.newick trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
		mv $treefile trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
	fi
fi

#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}MRL
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}MRL
fi

#Make MRP matrix
echo -e "Preparing MRL matrix...\n"
#java -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRPmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex NEXUS
#Make MRL matrix
if [[ $location == "2" ]]; then
	java -d64 -server -XX:MaxHeapSize=4g -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip PHYLIP -randomize
else
	java -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip PHYLIP -randomize
fi

#Make bootstrap ML trees using RAxML
if [[ $raxmlboot == "standard" ]]; then
	echo -e "Computing $bsrep standard bootstrap replicates using RAxML...\n"
	if [[ $location == "1" ]]; then
		raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n bestML -m BINCAT -p 12345
		raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -b 12345 -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n boot -m BINCAT -p 12345 -N $bsrep
		raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f b -t RAxML_bestTree.bestML -z RAxML_bootstrap.boot -n MRLresult -m BINCAT -p 12345
		mv RAxML_bootstrap.boot RAxML_bootstrap.MRLresult
		mv RAxML_bestTree.bestML RAxML_bestTree.MRLresult
		cat RAxML_info.bestML RAxML_info.boot RAxML_info.MRLresult > tmp && mv tmp RAxML_info.MRLresult
	elif [[ $location == "2" ]]; then
		raxmlHPC-PTHREADS-SSE3 -T $NSLOTS -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n bestML -m BINCAT -p 12345
		raxmlHPC-PTHREADS-SSE3 -T $NSLOTS -b 12345 -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n boot -m BINCAT -p 12345 -N $bsrep
		raxmlHPC-PTHREADS-SSE3 -T $NSLOTS -f b -t RAxML_bestTree.bestML -z RAxML_bootstrap.boot -n MRLresult -m BINCAT -p 12345
		mv RAxML_bootstrap.boot RAxML_bootstrap.MRLresult
		mv RAxML_bestTree.bestML RAxML_bestTree.MRLresult
		cat RAxML_info.bestML RAxML_info.boot RAxML_info.MRLresult > tmp && mv tmp RAxML_info.MRLresult
	elif [[ $numbcores == "1" ]]; then
		$raxmlseq -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n bestML -m BINCAT -p 12345
		$raxmlseq -b 12345 -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n boot -m BINCAT -p 12345 -N $bsrep
		$raxmlseq -f b -t RAxML_bestTree.bestML -z RAxML_bootstrap.boot -n MRLresult -m BINCAT -p 12345
		mv RAxML_bootstrap.boot RAxML_bootstrap.MRLresult
		mv RAxML_bestTree.bestML RAxML_bestTree.MRLresult
		cat RAxML_info.bestML RAxML_info.boot RAxML_info.MRLresult > tmp && mv tmp RAxML_info.MRLresult
	else
		$raxmlpthreads -T $numbcores -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n bestML -m BINCAT -p 12345
		$raxmlpthreads -T $numbcores -b 12345 -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n boot -m BINCAT -p 12345 -N $bsrep
		$raxmlpthreads -T $numbcores -f b -t RAxML_bestTree.bestML -z RAxML_bootstrap.boot -n MRLresult -m BINCAT -p 12345
		mv RAxML_bootstrap.boot RAxML_bootstrap.MRLresult
		mv RAxML_bestTree.bestML RAxML_bestTree.MRLresult
		cat RAxML_info.bestML RAxML_info.boot RAxML_info.MRLresult > tmp && mv tmp RAxML_info.MRLresult
	fi
else #rapid bootstrap
	echo -e "Computing $bsrep rapid bootstrap replicates using RAxML...\n"
	if [[ $location == "1" ]]; then
		raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N $bsrep
	elif [[ $location == "2" ]]; then
		raxmlHPC-PTHREADS-SSE3 -T $NSLOTS -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N $bsrep
	elif [[ $numbcores == "1" ]]; then
		$raxmlseq -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N $bsrep
	else
		$raxmlpthreads -T $numbcores -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N $bsrep
	fi
fi

#Modify labels in RAxML bipartitions (XX and YY to ' ')
sed -i.bak 's/XX/-/g' RAxML_bipartitions.MRLresult
sed -i.bak 's/YY/_/g' RAxML_bipartitions.MRLresult

#(Re)root/rename a final MRL species tree with $OUTGROUP
if [ -n "$OUTGROUP" ]; then
	nw_reroot -s RAxML_bipartitions.MRLresult $OUTGROUP > MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
else
	cp RAxML_bipartitions.MRLresult MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
fi

#Modify labels in RAxML bipartitions (XX and YY to ' ')
sed -i.bak 's/-/ /g' MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
sed -i.bak2 's/_/ /g' MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
sed -i.bak3 's/e 0/e-0/g' MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre

#Delete all *.bak files
rm *.bak*

#Rename/delete files
mv RAxML_info.MRLresult RAxML_MRL_info.log
mv RAxML_bootstrap.MRLresult MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre
rm RAxML_bipartitionsBranchLabels.MRLresult
rm RAxML_bipartitions.MRLresult
rm RAxML_bestTree.MRLresult

#Rename MRL trees
if [[ $requisite =~ "yes" ]]; then
	if [[ ! $collapse -eq "0" ]]; then
		mrltree=MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite_collapsed${collapse}
	else
		mrltree=MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite
	fi
else
	if [[ ! $collapse -eq "0" ]]; then
		mrltree=MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}_collapsed${collapse}
	else
		mrltree=MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}
	fi
fi

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp RAxML_MRL_info.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}MRL/${mrltree}_info.log
	cp MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}MRL/${mrltree}_allbootstraptrees.tre
	cp MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}MRL/${mrltree}.tre
else
	cp RAxML_MRL_info.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}MRL/${mrltree}_info.log
	cp MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}_allbootstraptrees.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}MRL/${mrltree}_allbootstraptrees.tre
	cp MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}MRL/${mrltree}.tre
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08c
fi

echo -e "HybPhyloMaker8c finished...\n"
