#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
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
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                    Script 06a - RAxML gene tree building                     *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Compute ML gene trees using RAxML for selected genes
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
	raxmlseq=raxmlHPC
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
	#Remove possible previously generated file for jobs
	rm -f ../submitRAxMLjobs.sh
	#Create new 'submitRAxMLjobs.sh' and make it executable
	touch ../submitRAxMLjobs.sh
	chmod +x ../submitRAxMLjobs.sh
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

#Check compatible setting (corrected=no is incompatible with genepart=codon)
if [[ $genetreepart == "codon" ]] && [[ $corrected == "no" ]]; then
	echo "You have incompatible settings (partitioning by codon [genetreepart=codon] is not allowed with uncorrected data [corrected=no]."
	echo "Change the settings before running the script again..."
	rm -d ../workdir06a/ 2>/dev/null
	exit 3
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -d "$path/$alnpathselected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}" ]; then
	if [ "$(ls -A $path/$alnpathselected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT})" ]; then
		if [ -f "$path/$alnpathselected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
			echo -e "OK\n"
		else
			echo -e "'$path/$alnpathselected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
			rm -d ../workdir06a/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/$alnpathselected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is empty. Exiting...\n"
		rm -d ../workdir06a/ 2>/dev/null
		exit 3
	fi
else
	echo -e "'$path/$alnpathselected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is missing. Exiting...\n"
	rm -d ../workdir06a/ 2>/dev/null
	exit 3
fi

#Test if folder for results exits
if [ -d "$path/$treepath${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML" ]; then
	echo -e "Directory '$path/$treepath${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML' already exists. Delete it or rename before running this script again. Exiting...\n"
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

#Write log
logname=HPM6a
echo -e "HybPhyloMaker6a: RAxML gene tree building" > ${logname}.log
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
for set in data selection cp corrected MISSINGPERCENT SPECIESPRESENCE genetreepart model raxmlboot bsrep bootstop numbcores raxmlperjob; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Add necessary scripts and files
cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}
mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
#Create new 'submitRAxMLjobs.sh' and make it executable
touch $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/submitRAxMLjobs.sh
chmod +x $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/submitRAxMLjobs.sh

if [[ $location == "1" || $location == "2" ]]; then
	echo -e "\nGenerating multiple jobs with $raxmlperjob alignments per job..."
	#Run on cluster: generate many jobs, each calculation only $raxmlperjob trees
	#Divide selected_genes$CUT.txt into files by $raxmlperjob
	split --lines=$raxmlperjob selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.
	rm selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
	cp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.* $path/${alnpathselected}${MISSINGPERCENT}
	for group in $(ls selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.*)
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
			echo 'qsub '"${group}"'.sh' >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML/submitRAxMLjobs.sh
		else
			cp ${group}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
			cp ${group}.sh ..
			echo 'qsub '"${group}"'.sh' >> ../submitRAxMLjobs.sh
		fi
	done
else
	#Run locally, trees are generated serially one by one
	if [[ $bootstop == "no" ]]; then
		echo -e "\nGenerating RAxML trees with $bsrep $raxmlboot bootstrap replicates using $model model...\n"
	else
		echo -e "\nGenerating RAxML trees with $raxmlboot bootstrap replicates using bootstop approach and $model model...\n"
	fi
	# Copy and modify selected FASTA files
	echo -e "Modifying selected FASTA files...\n"
	for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
	do
		cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .
		if [[ $genetreepart == "exon" ]]; then
			cp $path/${alnpath}/${i}.part .
		elif [[ $genetreepart == "codon" ]]; then
			cp $path/${alnpath}/${i}.codonpart.file .
			mv ${i}.codonpart.file ${i}.part
		fi
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
		#1.Check if there are completely undetermined columns in alignment (RAxML -y will produced .reduced alignment and partition files)
		#  Compute parsimony tree only and produce reduced alignment and appropriate reduced partition file
		if [[ $genetreepart == "no" ]]; then
			$raxmlseq -y -m $model -p 12345 -s $file.fas -n $file.check >> raxml.log
		else
			$raxmlseq -y -m $model -p 12345 -s $file.fas -q $filepart.part -n $file.check >> raxml.log
		fi
		#2.Test if reduced files were produced
		if [ -f $file.fas.reduced ]; then
			echo "Reduced alignment found...using it" >> raxml.log
			echo "Reduced alignment found...using it"
			#Put identical sequences back to the alignment (removed by RAxML when producing *.reduced dataset)
			#i.e., final alignemnt will have undetermined position removed but all samples present)
			grep "exactly identical" RAxML_info.${file}.check | grep "WARNING" | awk -F "Sequences |and |are " '{print $2 $3}' > recombine.txt
			nrlines=$(cat recombine.txt | wc -l )
			cat recombine.txt | while read -r a b
			do
				grep $a $file.fas.reduced > toadd.txt
				sed -i "s/$a/$b/" toadd.txt
				cat $file.fas.reduced toadd.txt > tmp && mv tmp $file.fas.reduced
			done
			rm recombine.txt toadd.txt 2>/dev/null
			#Correct number of samples in the final phylip file
			nrreduced=$(head -1 $file.fas.reduced | cut -d " " -f1)
			num=`expr $nrlines + $nrreduced`
			sed -i "1s/$nrreduced/$num/" $file.fas.reduced
			if [[ $genetreepart == "no" ]]; then
				#rm $file.fas
				mv $file.fas.reduced $file.fas
			else
				#rm $filepart.part
				mv $filepart.part.reduced $filepart.part
				#rm $file.fas
				mv $file.fas.reduced $file.fas
			fi
		else
			echo "Reduced alignment not found...using original alignment" >> raxml.log
			echo "Reduced alignment not found...using original alignment"
		fi
		#3.Run RAxML
		if [[ $raxmlboot == "standard" ]]; then
			if [[ $genetreepart == "no" ]]; then
				if [[ $bootstop == "no" ]]; then
					$raxmlpthreads -T $numbcores -s $file.fas -n $file.bestML -m $model -p 12345 >> raxml.log
					$raxmlpthreads -T $numbcores -b 12345 -s $file.fas -n $file.boot -m $model -p 12345 -N $bsrep >> raxml.log
					$raxmlpthreads -T $numbcores -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml.log
					mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result
					mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result
					cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result
				else
					$raxmlpthreads -T $numbcores -s $file.fas -n $file.bestML -m $model -p 12345 >> raxml.log
					$raxmlpthreads -T $numbcores -B 0.03 -b 12345 -s $file.fas -n $file.boot -m $model -p 12345 -N autoMRE >> raxml.log
					$raxmlpthreads -T $numbcores -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml.log
					bstop=$(grep "bootstrapped trees" RAxML_info.${file}.boot | awk '{ print $2 }')
					echo -e "$file\t$bstop" >> bootstop_summary.txt
					mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result
					mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result
					cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result
				fi
			else
				if [[ $bootstop == "no" ]]; then
					$raxmlpthreads -T $numbcores -s $file.fas -q $filepart.part -n $file.bestML -m $model -p 12345 >> raxml.log
					$raxmlpthreads -T $numbcores -b 12345 -s $file.fas -q $filepart.part -n $file.boot -m $model -p 12345 -N $bsrep >> raxml.log
					$raxmlpthreads -T $numbcores -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml.log
					mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result
					mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result
					cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result
				else
					$raxmlpthreads -T $numbcores -s $file.fas -q $filepart.part -n $file.bestML -m $model -p 12345 >> raxml.log
					$raxmlpthreads -T $numbcores -B 0.03 -b 12345 -s $file.fas -q $filepart.part -n $file.boot -m $model -p 12345 -N autoMRE >> raxml.log
					$raxmlpthreads -T $numbcores -f b -t RAxML_bestTree.${file}.bestML -z RAxML_bootstrap.${file}.boot -n $file.result -m $model -p 12345 >> raxml.log
					bstop=$(grep "bootstrapped trees" RAxML_info.${file}.boot | awk '{ print $2 }')
					echo -e "$file\t$bstop" >> bootstop_summary.txt
					mv RAxML_bootstrap.${file}.boot RAxML_bootstrap.${file}.result
					mv RAxML_bestTree.${file}.bestML RAxML_bestTree.${file}.result
					cat RAxML_info.${file}.bestML RAxML_info.${file}.boot RAxML_info.${file}.result > tmp && mv tmp RAxML_info.${file}.result
				fi
			fi
		else
			if [[ $genetreepart == "no" ]]; then
				if [[ $bootstop == "no" ]]; then
					$raxmlpthreads -T $numbcores -f a -s $file.fas -n $file.result -m $model -p 1234 -x 1234 -N $bsrep >> raxml.log
				else
					$raxmlpthreads -T $numbcores -B 0.03 -f a -s $file.fas -n $file.result -m $model -p 12345 -x 12345 -N autoMRE >> raxml.log
					bstop=$(grep "bootstrapped trees" RAxML_info.${file}.result | awk '{ print $2 }')
					echo -e "$file\t$bstop" >> bootstop_summary.txt
				fi
			else
				if [[ $bootstop == "no" ]]; then
					$raxmlpthreads -T $numbcores -f a -s $file.fas -q $filepart.part -n $file.result -m $model -p 1234 -x 1234 -N $bsrep >> raxml.log
				else
					$raxmlpthreads -T $numbcores -B 0.03 -f a -s $file.fas -q $filepart.part -n $file.result -m $model -p 12345 -x 12345 -N autoMRE >> raxml.log
					bstop=$(grep "bootstrapped trees" RAxML_info.${file}.result | awk '{ print $2 }')
					echo -e "$file\t$bstop" >> bootstop_summary.txt
				fi
			fi
		fi
		cp *${file}.result $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	done
	#Copy raxml.log to home
	cp raxml.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	if [[ $bootstop == "yes" ]]; then
		#Copy bootstop_summary.txt to home
		cp bootstop_summary.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML
	fi
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML

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
	echo -e "\nGo to homedir and run submitRAxMLjobs.sh...\n"
elif [[ $location == "1" ]]; then
	echo -e "\nGo to $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/RAxML and run submitRAxMLjobs.sh..."
	echo -e "This starts parallel computation of gene trees."
	echo -e "\nAfter all jobs finish run script HybPhyloMaker6a2 in order to calculate tree properties...\n"
elif [[ $location == "0" ]]; then
	echo -e "\nRun script HybPhyloMaker6a2 in order to calculate tree properties...\n"
fi
