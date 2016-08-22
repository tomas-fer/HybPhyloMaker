#!/bin/bash
#PBS -l walltime=2d
#PBS -l nodes=4:ppn=16:infiniband
#PBS -j oe
#PBS -l mem=4gb
#PBS -l scratch=32gb:shared
#PBS -N HybPhyloMaker7f_ExaMLconcatenated_tree
#PBS -m abe

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                     Script 07f - ExaML concatenated tree                     *
# *                                   v.1.1.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#IMPORTANT - submit to queue wagap.cerit-sc.cz using zuphux frontend (brno3-cerit)!!!
#unfortunately, may be it works from skirit as well (previously problem with RAxML within PartitionFinder)???

#Compute species tree from selected concatenated genes using by gene partitioning analysis in ExaML
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPhyloMaker4_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#or specify another input genes for concatenation below
#
#Calculate single best ExaML tree and than 100 bootstrap replicates (serially, could be slow...)
#Parallel version (by submitting multiple jobs) under development...

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
	module add fasttree-2.1.8
	module add python-3.4.1-intel
	module add examl-3.0.15
	module add raxml-8.2.8
	module add perl-5.20.1-gcc
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07f
	cd workdir07f
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

echo -e "\nHybPhyloMaker7f is running...\n"

#Check if there is already partition file (PartitionFinder was already run)
if [ -f $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/RAxMLpartitions.txt ]; then
	echo "Partition file found, skipping PartitionFinder run..."
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/RAxMLpartitions.txt .
	cp $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
else
	echo "Running concatenation and PartitionFinder first..."
	#Add necessary scripts and files
	cp $source/AMAS.py .
	#Copy list of genes
	if [[ $update =~ "yes" ]]; then
		cp $path/${alnpathselected}${type}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		mv selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
	else
		cp $path/${alnpathselected}${type}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
	fi
	
	# Make new dir for results
	if [[ $update =~ "yes" ]]; then
		mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
		mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	else
		mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	fi
	
	# Copy and modify selected FASTA files
	for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d"_" -f2); do
		#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
		cp $path/${alnpathselected}${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_${i}_modif${MISSINGPERCENT}.fas .
		#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
		sed -i.bak 's/(/_/g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
		sed -i.bak 's/)//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
		#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
		sed -i.bak 's/_contigs//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
		sed -i.bak 's/.fas//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	done
	
	#Prepare concatenated dataset and transform it to phylip format
	python3 AMAS.py concat -i *.fas -f fasta -d dna -u fasta -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta
	python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
	#Modify partition file: take 2nd and following parts of each line separated by "_" (remove 'p1_' etc. introduced by AMAS) | add ';' at the end of each line
	cat partitions.txt | cut -d"_" -f2- | awk '{ print $0 ";" }' | sed 's/-/ - /g' | sed 's/=/ = /g' > part.file
	#Copy concatenated file to home
	#Make new dir for results
	if [[ $update =~ "yes" ]]; then
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	else
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	fi
	
	# Prepare partition_finder.cfg
	echo -e "alignment = concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip;" >> partition_finder.cfg
	echo "branchlengths = linked;" >> partition_finder.cfg
	echo "models = GTR+G;" >> partition_finder.cfg
	echo "model_selection = AICc;" >> partition_finder.cfg
	echo "[data_blocks]" >> partition_finder.cfg
	cat partition_finder.cfg part.file > tmp && mv tmp partition_finder.cfg
	rm part.file
	echo "[schemes]" >> partition_finder.cfg
	echo "search = rclusterf;" >> partition_finder.cfg
	
	# Make a dir for PartitionFinder and move input files there
	mkdir PartitionFinder
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip PartitionFinder/
	mv partition_finder.cfg PartitionFinder/partition_finder.cfg
	
	#Get the newest version of PartitionFinder (2.0.0 from GitHub)
	#git clone https://github.com/brettc/partitionfinder
	#Copy working version of PartitionFinder 2.0.0.
	cp $source/partitionfinder-2.0.0-pre13.tar.gz .
	tar xfz partitionfinder-2.0.0-pre13.tar.gz
	mv partitionfinder-2.0.0-pre13 partitionfinder
	#Add python modules necessary for PartitionFinder
	module add python-2.7.10-gcc
	module add python27-modules-gcc
	module add hdf5-1.8.12-gcc
	# Run PartitionFinder
	#python partitionfinder/PartitionFinder.py PartitionFinder --raxml --rcluster-max 1000 --rcluster-percent 10
	python partitionfinder/PartitionFinder.py PartitionFinder --raxml
	#Prepare RAxML (ExaML) partition file from PartitionFinder results
	grep "DNA, Subset" PartitionFinder/analysis/best_scheme.txt > RAxMLpartitions.txt
	echo -e "\nPartitionFinder finished...\n"
	#Copy results to home
	if [[ $update =~ "yes" ]]; then
		cp -r PartitionFinder $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
		cp RAxMLpartitions.txt $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	else
		mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
		cp -r PartitionFinder $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
		cp RAxMLpartitions.txt $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	fi
fi

#Check if there are completely undetermined columns in alignment (RAxML -y will produced .reduced alignment and partition files)
#Compute parsimony tree only and produce reduced alignment and appropriate reduced partition file
echo -e "\nComputing parsimony tree...\n"
raxmlHPC -y -m GTRCAT -p 12345 -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -q RAxMLpartitions.txt -n check
#Test if reduced files were produced, copy them to home and rename to original names
if [ -f RAxMLpartitions.txt.reduced ]; then
	echo "Reduced alignment found...using it"
	if [[ $update =~ "yes" ]]; then
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.reduced $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
		cp RAxMLpartitions.txt.reduced $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	else
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.reduced $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
		cp RAxMLpartitions.txt.reduced $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	fi
	rm RAxMLpartitions.txt
	mv RAxMLpartitions.txt.reduced RAxMLpartitions.txt
	rm concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
	mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.reduced concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
else
	echo "Reduced alignment not found...using original alignment"
fi

#Run ExaML (single bestML only, no bootstrap)
#Set starting time
time=`date +%s`
#Prepare binary file for ExaML
echo -e "\nPreparing binary file for ExaML...\n"
parse-examl -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -m DNA -q RAxMLpartitions.txt -n concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.partitioned
#Prepare starting tree using RAxML
raxmlHPC -y -m GTRCAT -p 12345 -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n RandomStartingTree
#Run ExaML
echo -e "\nBuilding best ExaML tree...\n"
mpirun $examlbin -t RAxML_parsimonyTree.RandomStartingTree -m GAMMA -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.partitioned.binary -n examl.tre
#Measure time
timenow=`date +%s`
timespan=`echo "scale=2;($timenow-$time)/60" | bc`

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp ExaML* $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	echo -e "BestML done...in $timespan mins" >> $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML/ExaML_bootstrap_progress.txt
else
	cp ExaML* $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	echo -e "BestML done...in $timespan mins" >> $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/ExaML_bootstrap_progress.txt
fi

#Bootstrap ExaML
#generate 100 bootstrap replicates
echo -e "\nPreparing bootstrap replicates using RAxML...\n"
$raxmlseq -N 100 -b 12345 -f j -m GTRCAT -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -q RAxMLpartitions.txt -n REPS

for (( i=0; i<=99; i++ ))
do
	echo -e "\nBootstrap replicate $i...\n"
	#Set starting time
	time=`date +%s`
	#generate parsimony starting tree
	$raxmlseq -y -m GTRCAT -p 12345 -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i} -q RAxMLpartitions.txt.BS${i} -n T${i}
	#Prepare binary file for ExaML
	parse-examl -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i} -m DNA -q RAxMLpartitions.txt.BS${i} -n concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i}
	#Run ExaML
	mpirun $examlbin -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i}.binary -m GAMMA -t RAxML_parsimonyTree.T${i} -n BINF_${i}
	#Measure time
	timenow=`date +%s`
	timespan=`echo "scale=2;($timenow-$time)/60" | bc`
	if [[ $update =~ "yes" ]]; then
		echo -e "Bootstrap $i done...in $timespan mins" >> $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML/ExaML_bootstrap_progress.txt
		cp ExaML_result.BINF_${i} $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	else
		echo -e "Bootstrap $i done...in $timespan mins" >> $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/ExaML_bootstrap_progress.txt
		cp ExaML_result.BINF_${i} $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	fi
done

#Combine all bootstrap trees
echo "Combining bootstrap trees..."
cat ExaML_result.BINF* > ExaML_bootstrap.tre
#Map bootstrap values onto the single best tree
echo "Mapping bootstrap support values onto bestML tree"
$raxmlseq -f b -m GTRGAMMA -z ExaML_bootstrap.tre -t ExaML_result.examl.tre -n ExaML_Bootstrap

#Delete checkpoint files
rm *Check*

#Rename some files
mv ExaML_result.examl.tre ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre
mv ExaML_bootstrap.tre ExaML_bootstrap.trees
mv RAxML_bipartitions.ExaML_Bootstrap ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre

#Copy results to home
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML/ExaMLfiles
	cp *ExaML* $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML/ExaMLfiles
	cp ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	cp ExaML_bootstrap.trees $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
	cp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML
else
	mkdir $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/ExaMLfiles
	cp *ExaML* $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML/ExaMLfiles
	cp ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	cp ExaML_bootstrap.trees $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
	cp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre $path/${treepath}${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/concatenatedExaML
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07f
fi

echo -e "\nScript HybPhyloMaker 7f finished...\n"
