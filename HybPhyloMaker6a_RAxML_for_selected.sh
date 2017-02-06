#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -l scratch=1gb
#PBS -N HybPhyloMaker6a_RAxML_for_selected_parallel
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker6a_RAxML_for_selected_parallel
#$ -o HybPhyloMaker6a_RAxML_for_selected_parallel.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                    Script 06a - RAxML gene tree building                     *
# *                                   v.1.3.3                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute ML gene trees using RAxML (100 rapid bootstrap replicates) for selected genes
# Selection is based on maximum missing data per sample allowed ($MISSINGPERCENT) and minimum species percentage presence per assembly ($SPECIESPRESENCE)
# Edit those two values in settings.cfg
# Run first HybPhylomaker5_missingdataremoval.sh with the same combination of $MISSINGPERCENT and $SPECIESPRESENCE values
# If running locally, gene trees are produced serially (can be very SLOW with large alignment and lot of loci)
# If running on cluster, separate jobs are produced, number of jobs and number of loci per job is controlled by $raxmlperjob
# MetaCentrum runs all jobs automatically, on Hydra go to homedir and run submitRAxMLjobs.sh
# After all gene trees are generated run HybPhyloMaker6a2_RAxML_trees_summary.sh to calculate tree properties and plot graphs

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker6a is running on MetaCentrum..."
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
	echo -e "\nHybPhyloMaker6a is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir06a
	cd workdir06a
else
	echo -e "\nHybPhyloMaker6a is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir06a
	cd workdir06a
fi

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
else
	echo -e "Working with exons\n"
	type="exons"
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -d "$path/$type/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}" ]; then
	if [ "$(ls -A $path/$type/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT})" ]; then
		if [ -f "$path/$type/71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
			echo -e "OK\n"
		else
			echo -e "'$path/$type/71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
			rm -d ../workdir06a/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/$type/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is empty. Exiting...\n"
		rm -d ../workdir06a/ 2>/dev/null
		exit 3
	fi
else
	echo -e "'$path/$type/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is missing. Exiting...\n"
	rm -d ../workdir06a/ 2>/dev/null
	exit 3
fi

#Test if folder for results exits
if [ -d "$path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML" ]; then
	echo -e "Directory '$path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' already exists. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir06a/ 2>/dev/null
	exit 3
else
	if [[ ! $location == "1" ]]; then
		if [ "$(ls -A ../workdir06a)" ]; then
			echo -e "Directory 'workdir06a' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir06a/ 2>/dev/null
			exit 3
		fi
	fi
fi

#Add necessary scripts and files
cp $path/$type/71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
mkdir -p $path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}
mkdir $path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

if [[ $location == "1" || $location == "2" ]]; then
	echo -e "\nGenerating multiple jobs with $raxmlperjob alignments per job..."
	#Run on cluster: generate many jobs, each calculation only $raxmlperjob trees
	#Divide selected_genes$CUT.txt into files by $raxmlperjob
	split --lines=$raxmlperjob selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.
	rm selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
	cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.* $path/$type/71selected${MISSINGPERCENT}
	for group in $(ls selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.*)
	do
		echo '#!/bin/bash' >> ${group}.sh
		echo '#----------------MetaCentrum----------------' >> ${group}.sh
		echo '#PBS -l walltime=1d' >> ${group}.sh
		echo '#PBS -l nodes=1:ppn=4' >> ${group}.sh
		echo '#PBS -j oe' >> ${group}.sh
		echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> ${group}.sh
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
		echo '  mkdir workdir06_'"${group}" >> ${group}.sh
		echo '  cd workdir06_'"${group}" >> ${group}.sh
		echo 'fi' >> ${group}.sh
		echo 'path='"$path" >> ${group}.sh
		echo 'source='"$source" >> ${group}.sh
		echo 'MISSINGPERCENT='"$MISSINGPERCENT" >> ${group}.sh
		echo 'SPECIESPRESENCE='"$SPECIESPRESENCE" >> ${group}.sh
		echo 'type='"$type" >> ${group}.sh
		echo 'location='"$location" >> ${group}.sh
		echo 'genetreepart='"$genetreepart" >> ${group}.sh
		echo 'raxmlpthreads='"$raxmlpthreads" >> ${group}.sh
		echo 'cp '"$path"'/$type/71selected${MISSINGPERCENT}/'"$group"' .' >> ${group}.sh
		echo 'cp '"$source"'/catfasta2phyml.pl .' >> ${group}.sh
		echo 'for i in $(cat '"$group"')' >> ${group}.sh
		echo 'do' >> ${group}.sh
		echo '  cp '"$path"'/$type/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .' >> ${group}.sh
		echo '  cp $path/$type/70concatenated_exon_alignments/${i}.part .' >> ${group}.sh
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
		echo 'for file in $(cat FileForRAxML.txt); do' >> ${group}.sh
		echo '  #modify name for partition file (remove '_modif${MISSINGPERCENT}')' >> ${group}.sh
		echo '  filepart=$(sed "s/_modif${MISSINGPERCENT}//" <<< $file)' >> ${group}.sh
		echo '  #RAxML with 100 rapid bootstrap' >> ${group}.sh
		echo '  if [[ $location == "1" ]]; then' >> ${group}.sh
		echo '    if [[ $genetreepart == "yes" ]]; then' >> ${group}.sh
		echo '      raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f a -s $file.phylip -q $filepart.part -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100' >> ${group}.sh
		echo '    else' >> ${group}.sh
		echo '      raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f a -s $file.phylip -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100' >> ${group}.sh
		echo '    fi'  >> ${group}.sh
		echo '  elif [[ $location == "2" ]]; then' >> ${group}.sh
		echo '    if [[ $genetreepart == "yes" ]]; then' >> ${group}.sh
		echo '      $raxmlpthreads -T $NSLOTS -f a -s $file.phylip -q $filepart.part -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100' >> ${group}.sh
		echo '    else' >> ${group}.sh
		echo '      $raxmlpthreads -T $NSLOTS -f a -s $file.phylip -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100' >> ${group}.sh
		echo '    fi'  >> ${group}.sh
		echo '  fi' >> ${group}.sh
		echo '  cp *$file.result '"$path"'/$type/72trees'"${MISSINGPERCENT}_${SPECIESPRESENCE}"'/RAxML' >> ${group}.sh
		echo 'done' >> ${group}.sh
		echo '#Clean scratch/work directory' >> ${group}.sh
		echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${group}.sh
		echo '  #delete scratch' >> ${group}.sh
		echo '  rm -rf $SCRATCHDIR/*' >> ${group}.sh
		echo 'else' >> ${group}.sh
		echo '  cd ..' >> ${group}.sh
		echo '  rm -r workdir06_'"${group}" >> ${group}.sh
		echo 'fi' >> ${group}.sh
		
		chmod +x ${group}.sh
		if [[ $location == "1" ]]; then
			cp ${group}.sh $path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
			qsub ${group}.sh
		else
			cp ${group}.sh $path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
			cp ${group}.sh ..
			echo 'qsub '"${group}"'.sh' >> ../submitRAxMLjobs.sh
		fi
	done
else
	#Run locally, trees are generated serially one by one
	echo -e "\nGenerating RAxML trees with 100 rapid bootstrap replicates...\n"
	# Copy and modify selected FASTA files
	echo -e "Modifying selected FASTA files...\n"
	for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
	do
		cp $path/$type/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .
		cp $path/$type/70concatenated_exon_alignments/${i}.part .
		#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML)
		sed -i.bak 's/(/_/g' ${i}_modif${MISSINGPERCENT}.fas
		sed -i.bak 's/)//g' ${i}_modif${MISSINGPERCENT}.fas
		#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
		sed -i.bak 's/_contigs//g' ${i}_modif${MISSINGPERCENT}.fas
		sed -i.bak 's/.fas//g' ${i}_modif${MISSINGPERCENT}.fas
	done
	#Make a list of all fasta files
	ls *.fas | cut -d"." -f1 > FileForRAxMLTrees.txt
	echo -e "Generating RAxML trees..."
	numbertrees=$(cat FileForRAxMLTrees.txt | wc -l)
	calculating=0
	for file in $(cat FileForRAxMLTrees.txt); do
		calculating=$((calculating + 1))
		echo -e "\nProcessing file: ${file}" >> raxml.log
		echo -e "Processing file: ${file} ($calculating out of $numbertrees)"
		#modify name for partition file (remove '_modif${MISSINGPERCENT}'
		filepart=$(sed "s/_modif${MISSINGPERCENT}//" <<< $file)
		#run RAxML
		if [[ $genetreepart == "yes" ]]; then
			$raxmlpthreads -T $numbcores -f a -s $file.fas -q $filepart.part -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100 >> raxml.log
		else
			$raxmlpthreads -T $numbcores -f a -s $file.fas -q -n $file.result -m GTRCAT -p 1234 -x 1234 -N 100 >> raxml.log
		fi
		cp *$file.result $path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	done
#Copy raxml.log to home
cp raxml.log $path/$type/72trees${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir06a
fi

echo -e "\nHybPhyloMaker 6a finished..."
if [[ $location == "2" ]]; then
	echo -e "\nGo to homedir and run submitRAxMLjobs.sh..."
elif [[ $location == "1" ]]; then
	echo -e "\nAfter all jobs finish run script HybPhyloMaker6a2 in order to calculate tree properties..."
fi
echo -e "\nRun script HybPhyloMaker6a2 in order to calculate tree properties..."
