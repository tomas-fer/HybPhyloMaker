#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=1gb
#PBS -N HybPipe5a_RAxML_for_selected_parallel
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe5a_RAxML_for_selected_parallel
#$ -o HybPipe5a_RAxML_for_selected_parallel.log

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                    Script 05a - RAxML gene tree building                     *
# *                                   v.1.0.4                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute ML gene trees using RAxML (100 rapid bootstrap replicates) for selected genes
# Selection is based on maximum missing data allowed
# Edit CUT (cut off for missing data, i.e., maximum allowance for missing data) in settings.cfg
# Run first HybPipe4_missingdataremoval.sh with the same $CUT value

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
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
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir05a
	cd workdir05a
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

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	type="_cp"
else
	type=""
fi
#Add necessary scripts and files
cp $path/71selected${type}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
mkdir $path/72trees${type}${MISSINGPERCENT}_${SPECIESPRESENCE}
mkdir $path/72trees${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
#Divide selected_genes$CUT.txt into files by $raxmlperjob
split --lines=$raxmlperjob selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.
rm selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.* $path/71selected${type}${MISSINGPERCENT}
for group in $(ls selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.*)
do
	echo '#!/bin/bash' >> ${group}.sh
	echo '#----------------MetaCentrum----------------' >> ${group}.sh
	echo '#PBS -l walltime=1d' >> ${group}.sh
	echo '#PBS -l nodes=1:ppn=4' >> ${group}.sh
	echo '#PBS -j oe' >> ${group}.sh
	echo '#PBS -l mem=1gb' >> ${group}.sh
	echo '#PBS -N RAxML_for_'"${group}" >> ${group}.sh
	echo '#-------------------HYDRA-------------------' >> ${group}.sh
	echo '#$ -S /bin/bash' >> ${group}.sh
	echo '#$ -pe mthread 4' >> ${group}.sh
	echo '#$ -q sThC.q' >> ${group}.sh
	echo '#$ -l mres=4G,h_data=4G,h_vmem=4G' >> ${group}.sh
	echo '#$ -cwd' >> ${group}.sh
	echo '#$ -j y' >> ${group}.sh
	echo '#$ -N RAxML_for_'"${group}" >> ${group}.sh
	echo '#$ -o RAxML_for_'"${group}"'.log' >> ${group}.sh
	echo '#Complete path and set configuration for selected location' >> ${group}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${group}.sh
	echo '  . /packages/run/modules-2.0/init/bash' >> ${group}.sh
	echo '  #Add necessary modules' >> ${group}.sh
	echo '  module add raxml-8.2.4' >> ${group}.sh
	echo '  module add perl-5.10.1' >> ${group}.sh
	echo '  cd $SCRATCHDIR' >> ${group}.sh
	echo 'else' >> ${group}.sh
	echo '  #Add necessary modules' >> ${group}.sh
	echo '  module load bioinformatics/raxml/8.2.7' >> ${group}.sh
	echo '  mkdir workdir05_'"${group}" >> ${group}.sh
	echo '  cd workdir05_'"${group}" >> ${group}.sh
	echo 'fi' >> ${group}.sh
	echo 'path='"$path" >> ${group}.sh
	echo 'source='"$source" >> ${group}.sh
	echo 'MISSINGPERCENT='"$MISSINGPERCENT" >> ${group}.sh
	echo 'SPECIESPRESENCE='"$SPECIESPRESENCE" >> ${group}.sh
	echo 'type='"$type" >> ${group}.sh
	echo 'location='"$location" >> ${group}.sh
	echo 'cp '"$path"'/71selected${type}${MISSINGPERCENT}/'"$group"' .' >> ${group}.sh
	echo 'cp '"$source"'/catfasta2phyml.pl .' >> ${group}.sh
	echo 'for i in $(cat '"$group"')' >> ${group}.sh
	echo 'do' >> ${group}.sh
	echo '  cp '"$path"'/71selected${type}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .' >> ${group}.sh
	echo '  #Substitute '"'('"' by '"'_'"' and '"')'"' by nothing ('"'('"' and '"')'"' not allowed in RAxML)' >> ${group}.sh
	echo '  sed -i '"'s/(/_/g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
	echo '  sed -i '"'s/)//g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
	echo '  #Delete '"'_contigs'"' and '"'.fas'"' from labels (i.e., keep only genus-species_nr)' >> ${group}.sh
	echo '  sed -i '"'s/_contigs//g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
	echo '  sed -i '"'s/.fas//g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
	echo '  perl catfasta2phyml.pl ${i}_modif${MISSINGPERCENT}.fas > ${i}_modif${MISSINGPERCENT}.phylip' >> ${group}.sh
	echo 'done' >> ${group}.sh
	echo '#Make a list of all phylip files' >> ${group}.sh
	echo 'ls *.phylip | cut -d"." -f1 > FileForRAxML.txt' >> ${group}.sh
	echo 'for file in $(cat FileForRAxML.txt)' >> ${group}.sh
	echo 'do' >> ${group}.sh
	echo '  #RAxML with 100 rapid bootstrap' >> ${group}.sh
	echo '  if [[ $location == "1" ]]; then' >> ${group}.sh
	echo '    raxmlHPC-PTHREADS -T 4 -f a -s $file.phylip -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100' >> ${group}.sh
	echo '  else' >> ${group}.sh
	echo '    raxmlHPC-PTHREADS-SSE3 -T 4 -f a -s $file.phylip -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100' >> ${group}.sh
	echo '  fi' >> ${group}.sh
	echo '  cp *$file.result '"$path"'/72trees'"${type}${MISSINGPERCENT}_${SPECIESPRESENCE}"'/RAxML' >> ${group}.sh
	echo 'done' >> ${group}.sh
	echo '#Clean scratch/work directory' >> ${group}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${group}.sh
	echo '  #delete scratch' >> ${group}.sh
	echo '  rm -rf $SCRATCHDIR/*' >> ${group}.sh
	echo 'else' >> ${group}.sh
	echo '  cd ..' >> ${group}.sh
	echo '  rm -r workdir05_'"${group}" >> ${group}.sh
	echo 'fi' >> ${group}.sh
	
	chmod +x ${group}.sh
	if [[ $location == "1" ]]; then
		cp ${group}.sh $path/72trees${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
		qsub ${group}.sh
	else
		cp ${group}.sh $path/72trees${type}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
		cp ${group}.sh ..
		echo 'qsub '"${group}"'.sh' >> ../submitRAxMLjobs.sh
	fi
done

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir05a
fi
