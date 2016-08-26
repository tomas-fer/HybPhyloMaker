#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=4d
#PBS -l nodes=1:ppn=12
#PBS -j oe
#PBS -l mem=4gb
#PBS -l scratch=8gb
#PBS -N HybPhyloMaker7c_MRL
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 12
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker7c_MRL
#$ -o HybPhyloMaker7c_MRL.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                        Script 07c - MRL species tree                         *
# *                                   v.1.1.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree using MRL methods using RAxML from trees saved in single gene tree file (with *.newick suffix)
#Take trees from /concatenated_exon_alignments/selected${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick
#Run first
#(1) HybPhyloMaker4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker5a_RAxML_for_selected.sh or HybPhyloMaker5b_FastTree_for_selected.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(3) HybPhyloMaker6_roottrees.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#or specify another input trees below

#Complete path and set configuration for selected location
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
	module add jdk-1.6.0
	module add raxml-8.2.4
	module add newick-utils-1.6
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07c
	cd workdir07c
	#Add necessary modules
	module load java/1.7
	module load bioinformatics/raxml/8.2.7
	module load bioinformatics/newickutilities/0.0
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07c
	cd workdir07c
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
#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	type="_cp"
else
	type=""
fi

echo -e "\nHybPhyloMaker7c is running...\n"

#Add necessary programs and files
cp $source/mrp.jar .

#Copy genetree file
if [[ $update =~ "yes" ]]; then
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
else
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick .
fi

#Make dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/MRL
else
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/MRL
fi

#Make MRP matrix
echo -e "Preparing MRP matrix...\n"
#java -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRPmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.nex NEXUS
#Make MRL matrix
if [[ $location == "2" ]]; then
	java -d64 -server -XX:MaxHeapSize=10g -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip PHYLIP -randomize
else
	java -jar mrp.jar trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip PHYLIP -randomize
fi

#Make 100 fast bootstrap ML trees using RAxML
echo -e "Computing RAxML tree...\n"
if [[ $location == "1" ]]; then
	raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N 100
elif [[ $location == "2" ]]; then
	raxmlHPC-PTHREADS-SSE3 -T $NSLOTS -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N 100
elif [[ $numbcores == "1" ]]; then
	$raxmlseq -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N 100
else
	$raxmlpthreads -T $numbcores -f a -s MRLmatrix_${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n MRLresult -m BINCAT -p 1234 -x 1234 -N 100
fi

#Modify labels in RAxML bipartitions (XX and YY to ' ')
sed -i.bak 's/XX/-/g' RAxML_bipartitions.MRLresult
sed -i.bak 's/YY/_/g' RAxML_bipartitions.MRLresult

#(Re)root a final MRL species tree with $OUTGROUP
nw_reroot RAxML_bipartitions.MRLresult $OUTGROUP > MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre

#Modify labels in RAxML bipartitions (XX and YY to ' ')
sed -i.bak 's/-/ /g' MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
sed -i.bak 's/_/ /g' MRL_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre

#Delete all *.bak files
rm *.bak

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp *MR* $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/MRL
else
	cp *MR* $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/MRL
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07c
fi

echo -e "HybPhyloMaker7c finished...\n"
