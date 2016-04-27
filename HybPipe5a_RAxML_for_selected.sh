#!/bin/bash
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=1gb
#PBS -N HybPipe5a_RAxML_for_selected_parallel
#PBS -m abe

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                    Script 05a - RAxML gene tree building                     *
# *                                   v.1.0.1                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute ML gene trees using RAxML (100 rapid bootstrap replicates) for selected genes
# Selection is based on maximum missing data allowed
# Edit CUT (cut off for missing data, i.e., maximum allowance for missing data) in settings.cfg
# Run first HybPipe4_missingdataremoval.sh with the same $CUT value

#Complete path and set configuration for selected location
if [ ! $LOGNAME == "" ]; then
	echo "Metacentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	
else
	echo "Local..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir05a
	cd workdir05a
fi

#Add necessary scripts and files
cp $path/71selected${CUT}/selected_genes${CUT}.txt .
mkdir $path/72trees${CUT}
mkdir $path/72trees${CUT}/RAxML
#Divide selected_genes$CUT.txt into files by $raxmlperjob
split --lines=$raxmlperjob selected_genes$CUT.txt selected_genes$CUT.
rm selected_genes$CUT.txt
cp selected_genes$CUT.* $path/71selected${CUT}
for group in $(ls selected_genes$CUT.*)
do
	echo '#!/bin/bash' >> ${group}.sh
	echo '#PBS -l walltime=1d' >> ${group}.sh
	echo '#PBS -l nodes=1:ppn=4' >> ${group}.sh
	echo '#PBS -j oe' >> ${group}.sh
	echo '#PBS -l mem=1gb' >> ${group}.sh
	echo '#PBS -N RAxML_for_'"${group}" >> ${group}.sh
	echo '. /packages/run/modules-2.0/init/bash' >> ${group}.sh
	echo '#Add necessary modules' >> ${group}.sh
	echo 'module add raxml-8.2.4' >> ${group}.sh
	echo 'module add perl-5.10.1' >> ${group}.sh
	echo 'cd $SCRATCHDIR' >> ${group}.sh
	echo 'path='"$path" >> ${group}.sh
	echo 'source='"$source" >> ${group}.sh
	echo 'CUT='"$CUT" >> ${group}.sh
	echo 'cp '"$path"'/71selected${CUT}/'"$group"' .' >> ${group}.sh
	echo 'cp '"$source"'/catfasta2phyml.pl .' >> ${group}.sh
	echo 'for i in $(cat '"$group"')' >> ${group}.sh
	echo 'do' >> ${group}.sh
	echo '  cp '"$path"'/70concatenated_exon_alignments/$i.fasta .' >> ${group}.sh
	echo '  #Substitute '"'('"' by '"'_'"' and '"')'"' by nothing ('"'('"' and '"')'"' not allowed in RAxML)' >> ${group}.sh
	echo '  sed -i '"'s/(/_/g'"' $i.fasta' >> ${group}.sh
	echo '  sed -i '"'s/)//g'"' $i.fasta' >> ${group}.sh
	echo '  #Delete '"'_contigs'"' and '"'.fas'"' from labels (i.e., keep only genus-species_nr)' >> ${group}.sh
	echo '  sed -i '"'s/_contigs//g'"' $i.fasta' >> ${group}.sh
	echo '  sed -i '"'s/.fas//g'"' $i.fasta' >> ${group}.sh
	echo '  perl catfasta2phyml.pl $i.fasta > $i.phylip' >> ${group}.sh
	echo 'done' >> ${group}.sh
	echo '#Make a list of all phylip files' >> ${group}.sh
	echo 'ls *.phylip | cut -d"." -f1 > FileForRAxML.txt' >> ${group}.sh
	echo 'for file in $(cat FileForRAxML.txt)' >> ${group}.sh
	echo 'do' >> ${group}.sh
	echo '  #RAxML with 100 rapid bootstrap' >> ${group}.sh
	echo '  raxmlHPC-PTHREADS -T 4 -f a -s $file.phylip -n $file.result -m GTRCAT -p 1234 -x 1234 -# 100' >> ${group}.sh
	echo '  cp *$file.result '"$path"'/72trees'"${CUT}"'/RAxML' >> ${group}.sh
	echo 'done' >> ${group}.sh
	chmod +x ${group}.sh
	cp ${group}.sh $path/72trees${CUT}/RAxML
	qsub ${group}.sh
done

#Clean scratch/work directory
if [ ! $LOGNAME == "" ]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir05a
fi
