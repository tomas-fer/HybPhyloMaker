#!/bin/bash
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -l mem=1gb
#PBS -N HybPipe7d_MP-EST
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                       Script 07d - MP-EST species tree                       *
# *                                   v.1.0.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Compute species tree using MP-EST from rooted trees saved in single gene tree file (with *.newick suffix)
#Take trees from /concatenated_exon_alignments/selected${CUT}RAxML/species_trees
#Run first
#(1) HybPipe4_missingdataremoval.sh with the same $CUT value
#(2) HybPipe5a_RAxML_for_selected.sh or HybPipe5b_FastTree_for_selected.sh with the same $CUT value
#(3) HybPipe6_roottrees.sh with the same $CUT value
#or specify another input trees below

#Complete path and set configuration for selected location
if [ ! $LOGNAME == "" ]; then
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
	module add mpest-1.5
	module add R-3.1.0
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir07d
	cd workdir07d
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

#Copy genetree file and R script for resolving polytomies
cp $path/${treepath}${CUT}/${tree}/species_trees/trees${CUT}_rooted_withoutBS.newick .
cp $source/resolvepolytomies.R .
#Resolve polytomies in gene tree file (polytomies not allowed by MP-EST)
name=trees${CUT}_rooted_withoutBS.newick
R --slave -f resolvepolytomies.R $name
#Make new dir for results
mkdir $path/${treepath}${CUT}/${tree}/species_trees/MP-EST
#Copy tree with randomly resolved polytomies to home
cp trees${CUT}_rooted_withoutBS.newick.bifurcating.newick $path/${treepath}${CUT}/${tree}/species_trees/MP-EST
#Make control file for MP-EST 1.5 (see manual for details)
#Works for trees without bootstrap support
rm -f control
#Print name of the gene tree file (the only file with *.newick suffix in the folder)
echo `ls *.bifurcating*` >> control
#Print '0' meaning 'do not calculate triple distance among trees'
echo "0" >> control
#Print '-1' meaning 'random seed'
echo "-1" >> control
#Get number of trees in a gene tree file (delete possible empty lines and count lines)
nrtrees=`cat *.bifurcating* | sed '/^$/d' | wc -l`
#Get number of species in first gene tree file (delete '(' and ')', replace ',' by newlines, delete ';' and empty lines
nrtaxa=`cat *.bifurcating* | head -n 1 | sed 's/[(|)]//g' | tr ',' '\n' | sed 's/;//g' | sed '/^$/d' | cut -d":" -f1 | wc -l`
#Print number of trees and number of taxa
echo "$nrtrees $nrtaxa" >> control
#Taxon names to the variable $taxa
taxa=`cat *.bifurcating* | head -n 1 | sed 's/[(|)]//g' | tr ',' '\n' | sed 's/;//g' | sed '/^$/d' | cut -d":" -f1`
#Print taxon name, '1', taxon name
for i in $taxa
do
	echo -e "$i\t1\t$i" >> control
done
echo "0" >> control

#Copy control file to home
cp control $path/${treepath}${CUT}/${tree}/species_trees/MP-EST

#Run MP-EST
mpest control

#Add translate block
#Test if translate block exists (if there is a word 'translate' $test will be 1, otherwise 0)
test=`grep -i 'translate' *.tre | wc -l`
if [ "$test" -eq 0 ]; then
	echo "Adding translate block from control"
	#Add first four lines from the resulting tree file
	cat *.tre | head -n 4 > first
	#Write a word "echo"
	echo "  translate" > second
	#Add a numbered list of species from 'control': delete first 4 lines using sed, take first column using cut,
	#awk print number before and ',' after that column, replace last comma by ';'
	#perl command explanation:
	#, - matches a comma
	#(?!.*,) - negative lookahead asserts that there wouldn't be a comma after that matched, so matching the last comma
	#s - DOTALL modifier makes dot to match even newline characters also
	cat control | sed '1,4d' | sed '$d' | cut -f1 | awk '{print "\t" NR " " $0 ","}' | perl -00pe 's/,(?!.*,)/;/s' > third
	#Add resulting tree, delete first four lines, delete last two lines (i.e., retain only lines with trees)
	cat *.tre | sed '1,4d' | head --lines=-2 > fourth
	#Add two last lines from resulting tree file
	cat *.tre | tail -n 2 > fifth
	#Add lines containing 'tree mpest', i.e. species tree
	grep 'tree mpest' *.tre > mpesttree
	#Combine all parts together and produce again NEXUS tree files
	#all trees
	cat first second third fourth fifth > allMP-EST_trees_deleted_percentage${CUT}.tre
	#only mpest species tree
	cat first second third mpesttree fifth > MP-EST_species_tree_deleted_percentage${CUT}.tre
	#rm first second third fourth fifth mpesttree
else
	echo "Translate block already present"
	cat *.tre > allMP-EST_trees_${CUT}.tre
	#Delete all lines with word 'tree' followed by space and any number (i.e., retain only 'tree mpest' tree within NEXUS file)
	sed '/tree [0-9*]/d' allMP-EST_trees_${CUT}.tre > MP-EST_species_tree_${CUT}.tre
fi

#Modify names in species tree
sed -i 's/XX/-/g' *MP*.tre
sed -i 's/YY/_/g' *MP*.tre

#Copy results to home
cp *.tre $path/${treepath}${CUT}/${tree}/species_trees/MP-EST

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	#rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir07d
fi
