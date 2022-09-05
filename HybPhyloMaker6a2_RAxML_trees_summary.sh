#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:0:0
#PBS -l select=1:ncpus=4:mem=4gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker6a2_RAxML_trees_summary
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker6a2_RAxML_trees_summary
#$ -o HybPhyloMaker6a2_RAxML_trees_summary.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                  Script 06a2 - summary of RAxML gene trees                   *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute summary for already generated gene trees with RAxML
# Run first HybPhyloMaker5_missingdataremoval.sh and HybPhyloMaker6a_RAxML_for_selected.sh with the same settings

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker6a2 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add debian10-compat
	module add R-3.4.3-gcc
	#Set package library for R
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker6a2 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir06a2
	cd workdir06a2
	#Add necessary modules
	module load tools/R/3.4.1
else
	echo -e "\nHybPhyloMaker6a2 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir06a2
	cd workdir06a2
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
		echo -e "...with corrected reading frame\n"
	else
		mafftpath=$type/60mafft
		alnpath=$type/70concatenated_exon_alignments
		alnpathselected=$type/71selected
		treepath=$type/72trees
		echo -e "\n"
	fi
else
	if [[ $corrected =~ "yes" ]]; then
		mafftpath=$type/$selection/61mafft_corrected
		alnpath=$type/$selection/80concatenated_exon_alignments_corrected
		alnpathselected=$type/$selection/81selected_corrected
		treepath=$type/$selection/82trees_corrected
		echo -e "...with corrected reading frame...and for selection: $selection\n"
	else
		mafftpath=$type/$selection/60mafft
		alnpath=$type/$selection/70concatenated_exon_alignments
		alnpathselected=$type/$selection/71selected
		treepath=$type/$selection/72trees
		echo -e "...and for selection: $selection\n"
	fi
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
	if [ -d "$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}" ]; then
		if [ "$(ls -A $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT})" ]; then
			if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML" ]; then
				if [ "$(ls -A $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML)" ]; then
					echo -e "OK\n"
				else
					echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' is empty. Exiting...\n"
					rm -d ../workdir06a2/ 2>/dev/null
					exit 3
				fi
			else
				echo -e "'$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' is missing. Exiting...\n"
				rm -d ../workdir06a2/ 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is empty. Exiting...\n"
			rm -d ../workdir06a2/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is missing. Exiting...\n"
		rm -d ../workdir06a2/ 2>/dev/null
		exit 3
	fi
else
	echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
	rm -d ../workdir06a2/ 2>/dev/null
	exit 3
fi

#Test if folder for results exits
if [ -f "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/gene_properties.txt" ]; then
	echo -e "File '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/gene_properties.txt' already exists. You are probably going to overwrite previous results. Delete this file or rename before running this script again. Exiting...\n"
	rm -d ../workdir06a2/ 2>/dev/null
	exit 3
else
	if [[ ! $location == "1" ]]; then
		if [ "$(ls -A ../workdir06a2)" ]; then
			echo -e "Directory 'workdir06a2' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir06a2/ 2>/dev/null
			exit 3
		fi
	fi
fi

#----------------Checking number of alignments and trees (and stop if not identical)----------------
echo -e "Checking number of alignments and RAxML trees..."
nralignments=$(wc -l < $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
nrtrees=$(find $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML -maxdepth 1 -name "*bipartitions.*Assembly*" -exec ls {} + | wc -l)
echo -e "Number of selected alignments: $nralignments"
echo -e "Number of RAxML trees: $nrtrees"
if [ "$nralignments" -eq "$nrtrees" ]; then
	echo -e "Number of alignments and trees is the same. Continue with summary calculations...\n"
else
	echo -e "Number of alignments and trees does not match. Exiting...\n"
	#output alignment names
	cat $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | sort | sed '/^$/d' > alns.txt
	cp alns.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	find $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML -maxdepth 1 -name "*bipartitions.*Assembly*" -exec ls {} + | cut -d'.' -f2 | sed "s/_modif${MISSINGPERCENT}//" | sort | sed '/^$/d' > trees.txt
	cp trees.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	diff --new-line-format="" --unchanged-line-format="" alns.txt trees.txt > missingtrees.txt
	cp missingtrees.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	if [[ $location == "1" || $location == "2" ]]; then
		echo -e "Preparing job files for creating missing trees (one tree per job)\n"
		#Create new 'missingRAxMLjobs.sh' and make it executable
		touch $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/missingRAxMLjobs.sh
		chmod +x $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/missingRAxMLjobs.sh
		split --lines=1 missingtrees.txt missingtrees.
		rm missingtrees.txt
		cp missingtrees.* $path/${alnpathselected}${MISSINGPERCENT}
		for group in $(ls missingtrees.*)
		do
			echo '#!/bin/bash' >> ${group}.sh
			echo '#----------------MetaCentrum----------------' >> ${group}.sh
			echo '#PBS -l walltime=24:0:0' >> ${group}.sh
			echo '#PBS -l select=1:ncpus=4:mem=1gb:scratch_local=1gb' >> ${group}.sh
			echo '#PBS -j oe' >> ${group}.sh
			echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> ${group}.sh
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
			echo '  module load bioinformatics/raxml/8.2.11' >> ${group}.sh
			echo '  mkdir workdir06_'"${group}" >> ${group}.sh
			echo '  cd workdir06_'"${group}" >> ${group}.sh
			echo 'fi' >> ${group}.sh
			echo 'path='"$path" >> ${group}.sh
			echo 'source='"$source" >> ${group}.sh
			echo 'MISSINGPERCENT='"$MISSINGPERCENT" >> ${group}.sh
			echo 'SPECIESPRESENCE='"$SPECIESPRESENCE" >> ${group}.sh
			echo 'type='"$type" >> ${group}.sh
			echo 'corrected='"$corrected" >> ${group}.sh
			echo 'selection='"$selection" >> ${group}.sh
			echo 'location='"$location" >> ${group}.sh
			echo 'genetreepart='"$genetreepart" >> ${group}.sh
			echo 'raxmlpthreads='"$raxmlpthreads" >> ${group}.sh
			echo 'raxmlseq='"$raxmlseq" >> ${group}.sh
			echo 'bsrep='"$bsrep" >> ${group}.sh
			echo 'model='"$model" >> ${group}.sh
			echo 'raxmlboot='"$raxmlboot" >> ${group}.sh
			echo 'bootstop='"$bootstop" >> ${group}.sh
			echo 'if [ -z "$selection" ]; then' >> ${group}.sh
			echo '  if [[ $corrected =~ "yes" ]]; then' >> ${group}.sh
			echo '    alnpath=$type/80concatenated_exon_alignments_corrected' >> ${group}.sh
			echo '    alnpathselected=$type/81selected_corrected' >> ${group}.sh
			echo '    treepath=$type/82trees_corrected' >> ${group}.sh
			echo '  else' >> ${group}.sh
			echo '    alnpath=$type/70concatenated_exon_alignments' >> ${group}.sh
			echo '    alnpathselected=$type/71selected' >> ${group}.sh
			echo '    treepath=$type/72trees' >> ${group}.sh
			echo '  fi' >> ${group}.sh
			echo 'else' >> ${group}.sh
			echo '  if [[ $corrected =~ "yes" ]]; then' >> ${group}.sh
			echo '    alnpath=$type/$selection/80concatenated_exon_alignments_corrected' >> ${group}.sh
			echo '    alnpathselected=$type/$selection/81selected_corrected' >> ${group}.sh
			echo '    treepath=$type/$selection/82trees_corrected' >> ${group}.sh
			echo '  else' >> ${group}.sh
			echo '    alnpath=$type/$selection/70concatenated_exon_alignments' >> ${group}.sh
			echo '    alnpathselected=$type/$selection/71selected' >> ${group}.sh
			echo '    treepath=$type/$selection/72trees' >> ${group}.sh
			echo '  fi' >> ${group}.sh
			echo 'fi' >> ${group}.sh
			echo 'cp '"$path"'/${alnpathselected}${MISSINGPERCENT}/'"$group"' .' >> ${group}.sh
			echo 'for i in $(cat '"$group"')' >> ${group}.sh
			echo 'do' >> ${group}.sh
			echo '  cp '"$path"'/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .' >> ${group}.sh
			echo '  if [[ $genetreepart == "exon" ]]; then' >> ${group}.sh
			echo '    cp $path/${alnpath}/${i}.part .' >> ${group}.sh
			echo '  elif [[ $genetreepart == "codon" ]]; then' >> ${group}.sh
			echo '    cp $path/${alnpath}/${i}.codonpart.file .' >> ${group}.sh
			echo '    mv ${i}.codonpart.file ${i}.part' >> ${group}.sh
			echo '  fi' >> ${group}.sh
			echo '  #Substitute '"'('"' by '"'_'"' and '"')'"' by nothing ('"'('"' and '"')'"' not allowed in RAxML)' >> ${group}.sh
			echo '  sed -i '"'s/(/_/g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
			echo '  sed -i '"'s/)//g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
			echo '  #Delete '"'_contigs'"' and '"'.fas'"' from labels (i.e., keep only genus-species_nr)' >> ${group}.sh
			echo '  sed -i '"'s/_contigs//g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
			echo '  sed -i '"'s/.fas//g'"' ${i}_modif${MISSINGPERCENT}.fas' >> ${group}.sh
			echo 'done' >> ${group}.sh
			echo '#Make a list of all *.fas files' >> ${group}.sh
			echo 'ls *.fas | cut -d"." -f1 > FileForRAxML.txt' >> ${group}.sh
			echo 'for file in $(cat FileForRAxML.txt); do' >> ${group}.sh
			echo '  #modify name for partition file (remove '_modif${MISSINGPERCENT}')' >> ${group}.sh
			echo '  filepart=$(sed "s/_modif${MISSINGPERCENT}//" <<< $file)' >> ${group}.sh
			echo '  echo "Analysing $filepart" >> raxml'"$group"'.log' >> ${group}.sh
			echo '  #RAxML with bootstrap' >> ${group}.sh
			echo '  #1.Check if there are completely undetermined columns in alignment (RAxML -y will produced .reduced alignment and partition files)' >> ${group}.sh
			echo '  #  Compute parsimony tree only and produce reduced alignment and appropriate reduced partition file' >> ${group}.sh
			echo '  if [[ $genetreepart == "no" ]]; then' >> ${group}.sh
			echo '    $raxmlseq -y -m $model -p 12345 -s $file.fas -n $file.check >> raxml'"$group"'.log' >> ${group}.sh
			echo '  else' >> ${group}.sh
			echo '    $raxmlseq -y -m $model -p 12345 -s $file.fas -q $filepart.part -n $file.check >> raxml'"$group"'.log' >> ${group}.sh
			echo '  fi' >> ${group}.sh
			echo '  #2.Test if reduced files were produced' >> ${group}.sh
			echo '  if [ -f $file.fas.reduced ]; then' >> ${group}.sh
			echo '    echo "Reduced alignment found...using it" >> raxml'"$group"'.log' >> ${group}.sh
			echo '    #Put identical sequences back to the alignment (removed by RAxML when producing *.reduced dataset)' >> ${group}.sh
			echo '    #i.e., final alignemnt will have undetermined position removed but all samples present)' >> ${group}.sh
			echo "    grep \"exactly identical\" RAxML_info."'${file}'".check | grep \"WARNING\" | awk -F \"Sequences |and |are \" '{print \$2 \$3}' > recombine.txt" >> ${group}.sh
			echo '    nrlines=$(cat recombine.txt | wc -l )' >> ${group}.sh
			echo '    cat recombine.txt | while read -r a b' >> ${group}.sh
			echo '    do' >> ${group}.sh
			echo '      grep $a $file.fas.reduced > toadd.txt' >> ${group}.sh
			echo '      sed -i "s/$a/$b/" toadd.txt' >> ${group}.sh
			echo '      cat $file.fas.reduced toadd.txt > tmp && mv tmp $file.fas.reduced' >> ${group}.sh
			echo '    done' >> ${group}.sh
			echo '    rm recombine.txt toadd.txt 2>/dev/null' >> ${group}.sh
			echo '    #Correct number of samples in the final phylip file' >> ${group}.sh
			echo '    nrreduced=$(head -1 $file.fas.reduced | cut -d " " -f1)' >> ${group}.sh
			echo '    num=`expr $nrlines + $nrreduced`' >> ${group}.sh
			echo '    sed -i "1s/$nrreduced/$num/" $file.fas.reduced' >> ${group}.sh
			echo '    if [[ $genetreepart == "no" ]]; then' >> ${group}.sh
			echo '      #rm $file.fas' >> ${group}.sh
			echo '      mv $file.fas.reduced $file.fas' >> ${group}.sh
			echo '    else' >> ${group}.sh
			echo '      #rm $filepart.part' >> ${group}.sh
			echo '      mv $filepart.part.reduced $filepart.part' >> ${group}.sh
			echo '      #rm $file.fas' >> ${group}.sh
			echo '      mv $file.fas.reduced $file.fas' >> ${group}.sh
			echo '    fi' >> ${group}.sh
			echo '  else' >> ${group}.sh
			echo '    echo "Reduced alignment not found...using original alignment" >> raxml'"$group"'.log' >> ${group}.sh
			echo '  fi' >> ${group}.sh
			echo '  #3.Run RAxML' >> ${group}.sh
			echo '  if [[ $raxmlboot == "standard" ]]; then' >> ${group}.sh
			echo '    if [[ $location == "1" ]]; then' >> ${group}.sh
			echo '      if [[ $genetreepart == "no" ]]; then' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -s $file.fas -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -b 12345 -s $file.fas -n $file.boot -m $model -p 12345 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result' >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result' >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -s $file.fas -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -B 0.03 -b 12345 -s $file.fas -n $file.boot -m $model -p 12345 -N autoMRE >> raxml'"$group"'.log' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.boot | awk '"'"'{ print $2 }'"'"')' >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt' >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result' >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result' >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        fi' >> ${group}.sh
			echo '      else' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -s $file.fas -q $filepart.part -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -b 12345 -s $file.fas -q $filepart.part -n $file.boot -m $model -p 12345 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result' >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result' >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        else'  >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -s $file.fas -q $filepart.part -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log'  >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -B 0.03 -b 12345 -s $file.fas -q $filepart.part -n $file.boot -m $model -p 12345 -N autoMRE >> raxml'"$group"'.log'  >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log'  >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.boot | awk '"'"'{ print $2 }'"'"')'  >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt'  >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result'  >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result'  >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      fi'  >> ${group}.sh
			echo '    elif [[ $location == "2" ]]; then' >> ${group}.sh
			echo '      if [[ $genetreepart == "no" ]]; then' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -s $file.fas -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -b 12345 -s $file.fas -n $file.boot -m $model -p 12345 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result' >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result' >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -s $file.fas -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -B 0.03 -b 12345 -s $file.fas -n $file.boot -m $model -p 12345 -N autoMRE >> raxml'"$group"'.log' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.boot | awk '"'"'{ print $2 }'"'"')' >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt' >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result' >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result' >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      else' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -s $file.fas -q $filepart.part -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -b 12345 -s $file.fas -q $filepart.part -n $file.boot -m $model -p 12345 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log' >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result' >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result' >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -s $file.fas -q $filepart.part -n $file.bestML -m $model -p 12345 >> raxml'"$group"'.log'  >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -B 0.03 -b 12345 -s $file.fas -q $filepart.part -n $file.boot -m $model -p 12345 -N autoMRE >> raxml'"$group"'.log'  >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml'"$group"'.log'  >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.boot | awk '"'"'{ print $2 }'"'"')'  >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt'  >> ${group}.sh
			echo '          mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result'  >> ${group}.sh
			echo '          mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result'  >> ${group}.sh
			echo '          cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      fi'  >> ${group}.sh
			echo '    fi' >> ${group}.sh
			echo '    cp *$file.result '"${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}"'/RAxML' >> ${group}.sh
			echo '    if [[ $bootstop == "yes" ]]; then' >> ${group}.sh
			echo '      cp bootstop_summary_'"$group"'.txt '"${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}"'/RAxML' >> ${group}.sh
			echo '    fi' >> ${group}.sh
			echo '  else' >> ${group}.sh
			echo '    if [[ $location == "1" ]]; then' >> ${group}.sh
			echo '      if [[ $genetreepart == "no" ]]; then' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f a -s $file.fas -n $file.result -m $model -p 1234 -x 1234 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -B 0.03 -f a -s $file.fas -n $file.result -m $model -p 12345 -x 12345 -N autoMRE >> raxml'"$group"'.log' >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.result | awk '"'"'{ print $2 }'"'"')' >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      else' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -f a -s $file.fas -q $filepart.part -n $file.result -m $model -p 1234 -x 1234 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          raxmlHPC-PTHREADS -T $TORQUE_RESC_TOTAL_PROCS -B 0.03 -f a -s $file.fas -q $filepart.part -n $file.result -m $model -p 12345 -x 12345 -N autoMRE >> raxml'"$group"'.log' >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.result | awk '"'"'{ print $2 }'"'"')' >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      fi'  >> ${group}.sh
			echo '    elif [[ $location == "2" ]]; then' >> ${group}.sh
			echo '      if [[ $genetreepart == "no" ]]; then' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -f a -s $file.fas -n $file.result -m $model -p 1234 -x 1234 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -B 0.03 -f a -s $file.fas -n $file.result -m $model -p 12345 -x 12345 -N autoMRE >> raxml'"$group"'.log' >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.result | awk '"'"'{ print $2 }'"'"')' >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      else' >> ${group}.sh
			echo '        if [[ $bootstop == "no" ]]; then' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -f a -s $file.fas -q $filepart.part -n $file.result -m $model -p 1234 -x 1234 -N $bsrep >> raxml'"$group"'.log' >> ${group}.sh
			echo '        else' >> ${group}.sh
			echo '          $raxmlpthreads -T $NSLOTS -B 0.03 -f a -s $file.fas -q $filepart.part -n $file.result -m $model -p 12345 -x 12345 -N autoMRE >> raxml'"$group"'.log' >> ${group}.sh
			echo '          bstop=$(grep "bootstrapped trees" RAxML_info.${file}.result | awk '"'"'{ print $2 }'"'"')' >> ${group}.sh
			echo '          echo -e "$file\t$bstop" >> bootstop_summary_'"$group"'.txt' >> ${group}.sh
			echo '        fi'  >> ${group}.sh
			echo '      fi'  >> ${group}.sh
			echo '    fi' >> ${group}.sh
			echo '    cp *$file.result '"${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}"'/RAxML' >> ${group}.sh
			echo '    if [[ $bootstop == "yes" ]]; then' >> ${group}.sh
			echo '      cp bootstop_summary_'"$group"'.txt '"${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}"'/RAxML' >> ${group}.sh
			echo '    fi' >> ${group}.sh
			echo '  fi' >> ${group}.sh
			echo 'done' >> ${group}.sh
			echo 'cp raxml'"$group"'.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' >> ${group}.sh
			echo 'cp bootstop_summary_'"$group"'.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' >> ${group}.sh
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
				cp ${group}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
				#qsub ${group}.sh
				echo 'qsub '"${group}"'.sh' >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/missingRAxMLjobs.sh
			else
				cp ${group}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
				cp ${group}.sh ..
				echo 'qsub '"${group}"'.sh' >> ../missingRAxMLjobs.sh
			fi
		done
		if [[ $location == "2" ]]; then
			echo -e "\nGo to homedir and run missingRAxMLjobs.sh..."
			exit 3
		elif [[ $location == "1" ]]; then
			echo -e "\nGo to $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML and run missingRAxMLjobs.sh..."
			echo -e "This starts parallel computation of missing gene trees."
			echo -e "\nAfter all jobs finish run script HybPhyloMaker6a2 again in order to calculate tree properties..."
			exit 3
		fi
	else
		echo "Missing trees are given in $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/missingtrees.txt"
		rm -d ../workdir06a2/ 2>/dev/null
		exit 3
	fi
fi

#----------------Make a summary table with statistical properties for trees using R----------------
#Copy script
cp $source/tree_props.r .
cp $source/treepropsPlot.r .
cp $source/LBscores.R .
mkdir trees
mkdir alignments
#Copy all fasta alignments to subfolder 'alignments'
cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt); do
	cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas alignments/
done
#Copy all RAxML tree files (*.tre) with bootstrap values to subfolder 'trees'
#cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/*bipartitions.*Assembly* trees/
find $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML -maxdepth 1 -name "*bipartitions.*Assembly*" -exec cp -t trees {} + #to avoid 'Argument list too long' error
#Rename RAxML trees
cd trees
for i in *; do
	mv "$i" `basename "$i" .result | cut -f2 -d "."`.tre
done
cd ..

#Set a log file for R outputs/error messages
touch R.log

#Run R script for tree properties calculation and for plotting histograms of resulting values
echo -e "\nCalculating tree properties...\n"
echo -e "Calculating tree properties using tree_props.r\n" >> R.log
R --slave -f tree_props.r >> R.log 2>&1
#Run R script for calculation of LB score
echo -e "Calculating and parsing LB score...\n"
echo -e "\nCalculating and parsing LB score using LBscores.R\n" >> R.log
R --slave -f LBscores.R >> R.log 2>&1
#Parse script output (LBscores.csv)
echo -e "Taxon\tSum\tNrTrees\tMean" > LBscoresPerTaxon.txt
echo -e "Locus\tLBscoreSD" > LBscoresSDPerLocus.txt
for i in $(cat LBscores.csv | sed 1d | cut -d"," -f2 | sort | uniq); do
	awk -F',' -v val=$i 'BEGIN { n=0; sum=0} { if ($2 == val) { sum+=$4; n+=1 } } END { print val "\t" sum "\t" n "\t" sum/n }' LBscores.csv >> LBscoresPerTaxon.txt
done
for i in $(cat LBscores.csv | sed 1d | cut -d"," -f5 | sort | uniq); do
	grep $i LBscores.csv | awk -F',' -v val=$i '{ if ($5 == val) result=$1 } END { print val "\t" result }' >> LBscoresSDPerLocus.txt
done
cp LBscores.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
cp LBscoresPerTaxon.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
cp LBscoresSDPerLocus.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
#Combine 'LBscoresSDPerLocus.txt' with 'tree_stats_table.csv'
awk '{ print $2 }' LBscoresSDPerLocus.txt > tmp && mv tmp LBscoresSDPerLocus.txt
paste tree_stats_table.csv LBscoresSDPerLocus.txt | tr "\t" "," > tmp && mv tmp tree_stats_table.csv
#Replace 'NaN' by '0' (otherwise following plotting in R will not work)
sed -i.bak 's/NaN/0/g' tree_stats_table.csv
echo -e "Plotting boxplots/histograms for tree properties...\n"
echo -e "\nPlotting boxplots/histograms for tree properties using treepropsPlot.r\n" >> R.log
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f treepropsPlot.r
	R --slave -f treepropsPlot.r >> R.log 2>&1
else
	R --slave -f treepropsPlot.r >> R.log 2>&1
fi

#Copy results to home
cp tree_stats_table.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
cp *.png $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

#----------------Combine tree summary table with alignment summary and print comparison plots----------------
#Copy script
cp $source/plotting_correlations.R .
echo -e "Combining alignment and tree properties...\n"
#Copy alignment summary
cp $path/${alnpathselected}${MISSINGPERCENT}/summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
#***Modify alignment summary***
#Remove '_Assembly' (now locus names start with a number)
sed -i.bak 's/Assembly_//g' summarySELECTED*.txt
#Take first line as a header
head -n1 summarySELECTED*.txt > head.txt
#Remove first line | sort
sed 1d summarySELECTED*.txt | sort > summarySELECTED_sorted.txt
#Combine header and sorted file
cat head.txt summarySELECTED_sorted.txt > tmp && mv tmp summarySELECTED_sorted.txt
#***Modify tree summary***
#Change ',' to TAB
sed -i.bak 's/,/\t/g' tree_stats_table.csv
#Take first line as a header
head -n1 tree_stats_table.csv > head2.txt
#Remove first line | sort
sed 1d tree_stats_table.csv | sort > tree_stats_table_sorted.csv
#Combine header and sorted file
cat head2.txt tree_stats_table_sorted.csv > tmp && mv tmp tree_stats_table_sorted.csv
#Combine both files
paste summarySELECTED_sorted.txt tree_stats_table_sorted.csv > combined.txt

#Rename colums
sed -i.bak 's/Alignment_length/Aln_length/' combined.txt
sed -i.bak 's/Missing_percent/Missing_perc/' combined.txt
sed -i.bak 's/Proportion_parsimony_informative/Prop_pars_inf/' combined.txt
sed -i.bak 's/MstatX_entropy/Aln_entropy/' combined.txt
sed -i.bak 's/Average_bootstrap/Bootstrap/' combined.txt
sed -i.bak 's/Average_branch_length/Branch_length/' combined.txt
sed -i.bak 's/Avg_p_dist/P_distance/' combined.txt
sed -i.bak 's/Slope/Satur_slope/' combined.txt
sed -i.bak 's/R_squared/Satur_R_sq/' combined.txt
sed -i.bak 's/LBscoreSD/LBscore_SD/' combined.txt

#Run comparison plots for RAxML trees
echo -e "Plotting gene properties correlations for RAxML trees...\n"
echo -e "\nPlotting gene properties correlations for RAxML trees using plotting_correlations.R\n" >> R.log
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f plotting_correlations.R
	R --slave -f plotting_correlations.R >> R.log 2>&1
else
	R --slave -f plotting_correlations.R >> R.log 2>&1
fi
cp genes_corrs.* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
rm genes_corrs.*
mv combined.txt gene_properties.txt
cp gene_properties.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

#Copy R.log to home
cp R.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

#Combine bootstopping logs
if [[ $bootstop == "yes" ]]; then
	cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/bootstop_summary_*.txt .
	cat bootstop_summary_*.txt > bootstop_summary.txt
	cp bootstop_summary.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir06a2
fi

echo -e "HybPhyloMaker 6a2 finished...\n"
