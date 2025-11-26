#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker8m_PhyloNet_prepare
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8m_PhyloNet
#$ -o HybPhyloMaker8m_PhyloNet.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *             Script 08m - PhyloNet using maximum pseudo-likelihood            *
# *                                   v.1.8.0f                                   *
# *                                  Tomas Fer                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute PhyloNet networks using maximum pseudo-likelihood method from selected gene trees
#Take gene trees specified in trees${MISSINGPERCENT}_${SPECIESPRESENCE}_rooted_withoutBS.newick in 'species_trees' folder
#Take Astral species tree

#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#(2) HybPhyloMaker6b_FastTree_for_selected.sh or HybPhyloMaker6a_RAxML_for_selected.sh to create gene trees
#(3) HybPhyloMaker7_roottrees.sh to create the file containing all gene trees with BS values removed

#Requires
#PhyloNet.jar in 'HybSeqSource'

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8m is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	#module add debian9-compat
	#set nr processors
	#cpu=$TORQUE_RESC_PROC
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8m is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08m
	cd workdir08m
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #adds NewickUtilities
	module load tools/R/3.4.1
	#module load bioinformatics/newickutilities/0.0
	#julia
	#cairosvg
	#cpu=$NSLOTS
	#Remove possible previously generated file for jobs
	rm -f ../submitPhyloNetjobs.sh
	#Create new 'submitPhyloNetjobs.sh' and make it executable
	touch ../submitPhyloNetjobs.sh
	chmod +x ../submitPhyloNetjobs.sh
else
	echo -e "\nHybPhyloMaker8m is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08m
	cd workdir08m
	#cpu=4
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
echo -ne "\n\nTesting if gene trees are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile}" ]; then
		echo -e "OK\n"
	else
		echo -e "no gene trees file called '${modif}${treefile}' found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/'. Run first HybPhyloMaker7_roottrees.sh. Exiting..."
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
else
	if [ -f "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile}" ]; then
		echo -e "OK\n"
	else
		echo -e "no gene trees file called '${modif}${treefile}' found in '${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/'. Run first HybPhyloMaker7_roottrees.sh. Exiting..."
		rm -d ../workdir08l 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08m 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08m 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08m)" ]; then
		echo -e "Directory 'workdir08m' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08m 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8m
echo -e "HybPhyloMaker8m: PhyloNet network" > ${logname}.log
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
for set in data MISSINGPERCENT SPECIESPRESENCE selection cp corrected update tree hstart hmax; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Copy gene trees (file with all gene trees)
echo -e "\nCopying gene trees..."
if [[ $update =~ "yes" ]]; then
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}${treefile} .
	fi
else
	if [ -z "$OUTGROUP" ]; then
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	else
		cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}${treefile} .
	fi
fi

#Modify labels in gene tree
sed -i.bak 's/XX/-/g' $treefile
sed -i.bak2 's/YY/_/g' $treefile
#Removing '_cpDNA' from names
sed -i.bak 's/_cpDNA//g' $treefile

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet
fi

#Rename gene tree file etc.
mv $treefile trees.nwk
data1=$(echo $data | sed 's:.*/::') #everything after last slash, for naming purposes
nrh=$(($hmax-$hstart+1)) #number of different reticulation schemes
#Count number of trees (and decrease by one as they are counted from '0')
nrtrees=$(wc -l < trees.nwk)
nrtrees1=$((nrtrees - 1))

#Prepare tree NEXUS file
#begin commands
echo -e "#NEXUS\n\nBEGIN TREES;" > trees.nex
#add trees from NEWICK file
nr=0
cat trees.nwk | while read -r a; do
	echo -e "Tree gt${nr} = ${a}" >> trees.nex
	nr=$((nr + 1))
done
#end command
echo -e "END;" >> trees.nex

#Prepare PhyloNet input files
cpu=1 #hardcoded to '1' here, but will be replaced later by the submission script to all available CPUs
echo -e "hstart = $hstart"
echo -e "hmax = $hmax"
for pn in $(seq $hstart $hmax); do
	echo -e "NEXUS for $pn hybridization"
	cp trees.nex PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus
	echo -e "\nBEGIN PHYLONET;" >> PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus
	echo -e "InferNetwork_MPL (gt0-gt${nrtrees1}) ${pn} -pl ${cpu} -x ${numruns} -di ${data1}_${pn}_reticulations.txt;" >> PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus
	echo -e "END;" >> PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus
	if [[ $update =~ "yes" ]]; then
		mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}
		cp PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}
	else
		mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}
		cp PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}
	fi
done

#Create new 'submitPhyloNetjobs.sh' and make it executable
if [[ $update =~ "yes" ]]; then
	touch $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/submitPhyloNetjobs.sh
	chmod +x $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/submitPhyloNetjobs.sh
else
	touch $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/submitPhyloNetjobs.sh
	chmod +x $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/submitPhyloNetjobs.sh
fi

#Create PhyloNet job files
if [[ $location == "1" || $location == "2" ]]; then
	echo -e "\nGenerating multiple jobs for HybPhyloMaker with max nr. hybridizations from $hstart to $hmax..."
	for pn in $(seq $hstart $hmax); do
		echo -e "Generating PhyloNet job file for $pn hybridizations"
		echo '#!/bin/bash' >> PhyloNet_${pn}.sh
		echo '#----------------MetaCentrum----------------' >> PhyloNet_${pn}.sh
		echo '#PBS -l walltime=24:0:0' >> PhyloNet_${pn}.sh
		echo '#PBS -l select=1:ncpus=20:mem=100gb:scratch_local=8gb' >> PhyloNet_${pn}.sh
		echo '#PBS -j oe' >> PhyloNet_${pn}.sh
		echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> PhyloNet_${pn}.sh
		echo '#PBS -N PhyloNet_for_'"${pn}" >> PhyloNet_${pn}.sh
		echo '#-------------------HYDRA-------------------' >> PhyloNet_${pn}.sh
		echo '#$ -S /bin/bash' >> PhyloNet_${pn}.sh
		echo '#$ -pe mthread 12' >> PhyloNet_${pn}.sh
		echo '#$ -q sThC.q' >> PhyloNet_${pn}.sh
		echo '#$ -l mres=200G,h_data=200G,h_vmem=200G' >> PhyloNet_${pn}.sh
		echo '#$ -cwd' >> PhyloNet_${pn}.sh
		echo '#$ -j y' >> PhyloNet_${pn}.sh
		echo '#$ -N PhyloNet_for_'"${pn}" >> PhyloNet_${pn}.sh
		echo '#$ -o PhyloNet_for_'"${pn}"'.log' >> PhyloNet_${pn}.sh
		echo '#Complete path and set configuration for selected location' >> PhyloNet_${pn}.sh
		echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> PhyloNet_${pn}.sh
		echo '  #Add necessary modules' >> PhyloNet_${pn}.sh
		echo '  module add jdk-8' >> PhyloNet_${pn}.sh
		echo '  #set nr processors and memory size' >> PhyloNet_${pn}.sh
		echo '  cpu=$TORQUE_RESC_PROC' >> PhyloNet_${pn}.sh
		echo '  cpuavail=$(($cpu-1))' >> PhyloNet_${pn}.sh
		echo '  mem=$PBS_RESC_MEM' >> PhyloNet_${pn}.sh
		echo '  cd $SCRATCHDIR' >> PhyloNet_${pn}.sh
		echo 'else' >> PhyloNet_${pn}.sh
		echo '  #Add necessary modules' >> PhyloNet_${pn}.sh
		echo '  module load java/1.7' >> PhyloNet_${pn}.sh
		echo '  mkdir workdir08m_'"${pn}" >> PhyloNet_${pn}.sh
		echo '  cd workdir08m_'"${pn}" >> PhyloNet_${pn}.sh
		echo 'fi' >> PhyloNet_${pn}.sh
		echo 'path='"$path" >> PhyloNet_${pn}.sh
		echo 'source='"$source" >> PhyloNet_${pn}.sh
		echo 'MISSINGPERCENT='"$MISSINGPERCENT" >> PhyloNet_${pn}.sh
		echo 'SPECIESPRESENCE='"$SPECIESPRESENCE" >> PhyloNet_${pn}.sh
		echo 'type='"$type" >> PhyloNet_${pn}.sh
		echo 'corrected='"$corrected" >> PhyloNet_${pn}.sh
		echo 'selection='"$selection" >> PhyloNet_${pn}.sh
		echo 'location='"$location" >> PhyloNet_${pn}.sh
		echo 'update='"$update" >> PhyloNet_${pn}.sh
		echo 'tree='"$tree" >> PhyloNet_${pn}.sh
		echo 'nrtrees='"$nrtrees" >> PhyloNet_${pn}.sh
		echo 'data='"$data" >> PhyloNet_${pn}.sh
		echo 'if [ -z "$selection" ]; then' >> PhyloNet_${pn}.sh
		echo '  data1=$(echo $data | sed '"'s:.*/::'"')' >> PhyloNet_${pn}.sh
		echo 'else' >> PhyloNet_${pn}.sh
		echo '  data1=$selection' >> PhyloNet_${pn}.sh
		echo 'fi' >> PhyloNet_${pn}.sh
		echo 'pn='"$pn" >> PhyloNet_${pn}.sh
		echo 'if [ -z "$selection" ]; then' >> PhyloNet_${pn}.sh
		echo '  if [[ $corrected =~ "yes" ]]; then' >> PhyloNet_${pn}.sh
		echo '    alnpath=$type/80concatenated_exon_alignments_corrected' >> PhyloNet_${pn}.sh
		echo '    alnpathselected=$type/81selected_corrected' >> PhyloNet_${pn}.sh
		echo '    treepath=$type/82trees_corrected' >> PhyloNet_${pn}.sh
		echo '  else' >> PhyloNet_${pn}.sh
		echo '    alnpath=$type/70concatenated_exon_alignments' >> PhyloNet_${pn}.sh
		echo '    alnpathselected=$type/71selected' >> PhyloNet_${pn}.sh
		echo '    treepath=$type/72trees' >> PhyloNet_${pn}.sh
		echo '  fi' >> PhyloNet_${pn}.sh
		echo 'else' >> PhyloNet_${pn}.sh
		echo '  if [[ $corrected =~ "yes" ]]; then' >> PhyloNet_${pn}.sh
		echo '    alnpath=$type/$selection/80concatenated_exon_alignments_corrected' >> PhyloNet_${pn}.sh
		echo '    alnpathselected=$type/$selection/81selected_corrected' >> PhyloNet_${pn}.sh
		echo '    treepath=$type/$selection/82trees_corrected' >> PhyloNet_${pn}.sh
		echo '  else' >> PhyloNet_${pn}.sh
		echo '    alnpath=$type/$selection/70concatenated_exon_alignments' >> PhyloNet_${pn}.sh
		echo '    alnpathselected=$type/$selection/71selected' >> PhyloNet_${pn}.sh
		echo '    treepath=$type/$selection/72trees' >> PhyloNet_${pn}.sh
		echo '  fi' >> PhyloNet_${pn}.sh
		echo 'fi' >> PhyloNet_${pn}.sh
		echo '#Copy PhyloNet jar file' >> PhyloNet_${pn}.sh
		echo 'wget https://phylogenomics.rice.edu/media/PhyloNet.jar' >> PhyloNet_${pn}.sh
		#echo 'wget https://github.com/NakhlehLab/PhyloNet/releases/latest/download/PhyloNet.jar' >> PhyloNet_${pn}.sh
		echo '#Copy PhyloNet NEXUS input file' >> PhyloNet_${pn}.sh
		echo 'if [[ $update =~ "yes" ]]; then' >> PhyloNet_${pn}.sh
		echo '  cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}/PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus .' >> PhyloNet_${pn}.sh
		echo 'else' >> PhyloNet_${pn}.sh
		echo '  cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}/PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus .' >> PhyloNet_${pn}.sh
		echo 'fi' >> PhyloNet_${pn}.sh
		echo '#Modify nr of CPUs to all available' >> PhyloNet_${pn}.sh
		echo 'sed -i "s/pl 1/pl ${cpuavail}/" PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus' >> PhyloNet_${pn}.sh
		echo '#Run PhyloNet' >> PhyloNet_${pn}.sh
		echo 'java -Xmx${mem} -jar PhyloNet.jar PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus > ${data1}_${pn}_reticulations.txt' >> PhyloNet_${pn}.sh
		echo 'sed -ne '"'"'/Inferred/,$ p'"'"' ${data1}_${pn}_reticulations.txt > ${data1}_${pn}_reticulations.networks' >> PhyloNet_${pn}.sh
		echo 'grep "probability" ${data1}_${pn}_reticulations.networks | cut -d'"' '"' -f4 > ${pn}_logliks.txt' >> PhyloNet_${pn}.sh
		echo 'grep -A1 '"\"Network\""' ${data1}_${pn}_reticulations.networks | grep '"\"(\""' | grep -on '"\"1\.0\" | cut -d':' -f1 | uniq -c | awk '{"'$1=$1'"};1' | cut -d' ' -f1 > "'${pn}'"_nrbranches.txt" >> PhyloNet_${pn}.sh
		echo 'for i in $(seq 1 5); do echo ${pn}; done >> ${pn}_nrhyb.txt' >> PhyloNet_${pn}.sh
		echo 'for i in $(seq 1 5); do echo ${nrtrees}; done >> ${pn}_nrgtrees.txt' >> PhyloNet_${pn}.sh
		echo 'paste ${pn}_nrhyb.txt ${pn}_nrbranches.txt ${pn}_nrgtrees.txt ${pn}_logliks.txt > ${pn}_table.txt' >> PhyloNet_${pn}.sh
		echo 'rm ${pn}_nrhyb.txt ${pn}_nrbranches.txt ${pn}_nrgtrees.txt ${pn}_logliks.txt' >> PhyloNet_${pn}.sh
		echo '#Copy results to home' >> PhyloNet_${pn}.sh
		echo 'if [[ $update =~ "yes" ]]; then' >> PhyloNet_${pn}.sh
		echo '  mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}' >> PhyloNet_${pn}.sh
		echo '  cp *.{networks,txt} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/${pn}' >> PhyloNet_${pn}.sh
		echo 'else' >> PhyloNet_${pn}.sh
		echo '  mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}' >> PhyloNet_${pn}.sh
		echo '  cp *.{networks,txt} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/${pn}' >> PhyloNet_${pn}.sh
		echo 'fi' >> PhyloNet_${pn}.sh
		echo '#Clean scratch/work directory' >> PhyloNet_${pn}.sh
		echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> PhyloNet_${pn}.sh
		echo '  #delete scratch' >> PhyloNet_${pn}.sh
		echo '  rm -rf $SCRATCHDIR/*' >> PhyloNet_${pn}.sh
		echo 'else' >> PhyloNet_${pn}.sh
		echo '  cd ..' >> PhyloNet_${pn}.sh
		echo '  rm -r workdir08m_'"${pn}" >> PhyloNet_${pn}.sh
		echo 'fi' >> PhyloNet_${pn}.sh
		chmod +x PhyloNet_${pn}.sh
		if [[ $location == "1" ]]; then
			if [[ $update =~ "yes" ]]; then
				cp PhyloNet_${pn}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet
				echo 'qsub PhyloNet_'"${pn}"'.sh' >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet/submitPhyloNetjobs.sh
			else
				cp PhyloNet_${pn}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet
				echo 'qsub PhyloNet_'"${pn}"'.sh' >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet/submitPhyloNetjobs.sh
			fi
		else
			if [[ $update =~ "yes" ]]; then
				cp PhyloNet_${pn}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet
				cp PhyloNet_${pn}.sh ..
				echo 'qsub PhyloNet_'"${pn}"'.sh' >> ../submitPhyloNetjobs.sh
			else
				cp PhyloNet_${pn}.sh $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet
				cp PhyloNet_${pn}.sh ..
				echo 'qsub PhyloNet_'"${pn}"'.sh' >> ../submitPhyloNetjobs.sh
			fi
		fi
	done
else
	echo -e "\nRunning PhyloNet locally, serially for nr. hybridizations from $hstart to $hmax. This might take a lot of time"
	for pn in $(seq $hstart $hmax); do
		echo -e "PhyloNet for nr. hybridization: ${pn}"
		java -jar PhyloNet.jar PhyloNet_${data1}_${nrtrees}genetrees_${pn}reticulations.nexus > ${data1}_${pn}_reticulations.txt
		sed -ne '/Inferred/,$ p' ${data1}_${pn}_reticulations.txt > ${data1}_${pn}_reticulations.networks
		grep "probability" ${data1}_${pn}_reticulations.networks | cut -d' ' -f4 > ${pn}_logliks.txt
		grep -A1 "Network" ${data1}_${pn}_reticulations.networks | grep "(" | grep -on "1\.0" | cut -d':' -f1 | uniq -c | awk '{$1=$1};1' | cut -d' ' -f1 > ${pn}_nrbranches.txt
		for i in $(seq 1 5); do echo ${pn}; done >> ${pn}_nrhyb.txt
		for i in $(seq 1 5); do echo ${nrtrees}; done >> ${pn}_nrgtrees.txt
		paste ${pn}_nrhyb.txt ${pn}_nrbranches.txt ${pn}_nrgtrees.txt ${pn}_logliks.txt > ${pn}_table.txt
		rm ${pn}_nrhyb.txt ${pn}_nrbranches.txt ${pn}_nrgtrees.txt ${pn}_logliks.txt
	done
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08m
fi

echo -e "\nHybPhyloMaker 8m finished..."
if [[ $location == "2" ]]; then
	echo -e "\nGo to homedir and run submitPhyloNetjobs.sh...\n"
	echo -e "This starts parallel computation of PhyloNet networks with max nr. hybridization from $hstart to $hmax."
	echo -e "\nAfter all jobs finish run script HybPhyloMaker8m2 in order to summarize PhyloNet results...\n"
elif [[ $location == "1" ]]; then
	if [[ $update =~ "yes" ]]; then
		echo -e "\nGo to $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/PhyloNet and run submitRAxMLjobs.sh..."
	else
		echo -e "\nGo to $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/PhyloNet and run submitRAxMLjobs.sh..."
	fi
	echo -e "This starts parallel computation of PhyloNet networks with max nr. hybridization from $hstart to $hmax."
	echo -e "\nAfter all jobs finish run script HybPhyloMaker8m2 in order to summarize PhyloNet results...\n"
elif [[ $location == "0" ]]; then
	echo -e "\nRun script HybPhyloMaker8m2 in order to summarize PhyloNet results...\n"
fi
