#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=4gb:scratch_local=32gb
#PBS -j oe
#PBS -N HybPhyloMaker8f_ExaMLconcatenated_tree
#PBS -m abe

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                     Script 08f - ExaML concatenated tree                     *
# *                                   v.1.8.0b                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#requires partitionfinder-2.1.1.tar.gz in HybSeqSource

#IMPORTANT - submit to queue wagap.cerit-sc.cz using zuphux frontend (brno3-cerit)!!!
#unfortunately, may be it works from skirit as well (previously problem with RAxML within PartitionFinder)???

#Compute species tree from selected concatenated genes using by gene partitioning analysis in ExaML
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values
#or specify another input genes for concatenation below
#
#Calculate single best ExaML tree and than 100 bootstrap replicates (serially, could be slow...)
#Parallel version (by submitting multiple jobs) under development...

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8f is running on MetaCentrum..."
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
	module add python-3.4.1-gcc
	module add examl-3.0.15
	module add raxml-8.2.8
	module add newick-utils-13042016
	module add perl-5.20.1-gcc
	#module add debian8-compat
else
	echo -e "\nHybPhyloMaker8f is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08f
	cd workdir08f
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

if [[ $requisite =~ "yes" ]]; then
	echo -e "...and only with trees with requisite taxa present\n"
else
	echo -e "\n"
fi

#Collapsed trees have no meaning in this concatenated analysis
collapse=0

#Settings for requisite selection
if [[ $requisite =~ "yes" ]]; then
	modif=with_requisite/
	modif1=with_requisite
	treefile=trees_rooted_with_requisite.newick
else
	modif=""
	modif1=""
	treefile=trees_rooted.newick
fi

#Check if there is already partition file (PartitionFinder was already run)
if [[ $update =~ "yes" ]]; then
	if [ -f $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/RAxMLpartitions.txt ]; then
		echo "Partition file found, skipping PartitionFinder run..."
		if [[ $requisite =~ "yes" ]]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/RAxMLpartitions.txt .
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.phylip .
			mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.phylip concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
			runpf=no
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/RAxMLpartitions.txt .
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
			runpf=no
		fi
	else
		runpf=yes
	fi
else
	if [ -f $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/RAxMLpartitions.txt ]; then
		echo "Partition file found, skipping PartitionFinder run..."
		if [[ $requisite =~ "yes" ]]; then
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/RAxMLpartitions.txt .
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.phylip .
			mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.phylip concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
			runpf=no
		else
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/RAxMLpartitions.txt .
			cp $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip .
			runpf=no
		fi
	else
		runpf=yes
	fi
fi


#Check necessary file
if [[ $requisite =~ "no" ]]; then
	echo -ne "Testing if input data are available..."
	if [[ $update =~ "yes" ]]; then
		if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
			#if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
			if [ 0 -lt $(find $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT} -maxdepth 1 -name "*ssembly_*_modif${MISSINGPERCENT}.fas" -exec ls {} + 2>/dev/null | wc -w) ]; then #to avoid 'Argument list too long' error
				echo -e "OK\n"
			else
				echo -e "no alignmenet files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
				rm -d ../workdir08f 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
			rm -d ../workdir08f 2>/dev/null
			exit 3
		fi
	else
		if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
			#if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
			if [ 0 -lt $(find $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT} -maxdepth 1 -name "*ssembly_*_modif${MISSINGPERCENT}.fas" -exec ls {} + 2>/dev/null | wc -w) ]; then #to avoid 'Argument list too long' error
				echo -e "OK\n"
			else
				echo -e "no alignmenet files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
				rm -d ../workdir08f 2>/dev/null
				exit 3
			fi
		else
			echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
			rm -d ../workdir08f 2>/dev/null
			exit 3
		fi
	fi
	
	if [[ $runpf =~ "yes" ]]; then
		#Test if folder for results exits
		if [[ $update =~ "yes" ]]; then
			if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML" ]; then
				echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenatedExaML' already exists. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir08f 2>/dev/null
				exit 3
			fi
		else
			if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenatedExaML" ]; then
				echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenatedExaML' already exists. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir08f 2>/dev/null
				exit 3
			fi
		fi
		if [[ ! $location == "1" ]]; then
			if [ "$(ls -A ../workdir08f)" ]; then
				echo -e "Directory 'workdir08f' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
				rm -d ../workdir08f 2>/dev/null
				exit 3
			fi
		fi
	fi
fi

#Write log
logname=HPM8f
echo -e "HybPhyloMaker8f: concatenated species tree using ExaML" > ${logname}.log
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
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree OUTGROUP requisite examlboot; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [[ $requisite =~ "yes" ]]; then
	echo -e "\nList of requisite samples" >> ${logname}.log
	echo $requisitetaxa | tr '|' '\n' >> ${logname}.log
fi
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

if [[ $runpf =~ "yes" ]]; then
	echo -e "Running concatenation and PartitionFinder first...\n"
	#Add necessary scripts and files
	cp $source/AMAS.py .
	#Copy list of genes
	if [[ $update =~ "yes" ]]; then
		cp $path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt .
		mv selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
	else
		cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
	fi
	
	# Make new dir for results
	if [[ $update =~ "yes" ]]; then
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	else
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	fi
	
	# Copy and modify selected FASTA files
	echo -e "Copying and modifying FASTA files...\n"
	for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d"_" -f2); do
		#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
		cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_${i}_modif${MISSINGPERCENT}.fas .
		if [[ $requisite =~ "yes" ]]; then
			#Delete an FASTA alignment if requisite taxa are not present in it (at least one of)
			if ! grep -E ${requisitetaxa} *ssembly_${i}_modif${MISSINGPERCENT}.fas >/dev/null; then
				rm *ssembly_${i}_modif${MISSINGPERCENT}.fas
			else
				#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
				sed -i.bak 's/(/_/g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
				sed -i.bak 's/)//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
				#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
				sed -i.bak 's/_contigs//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
				sed -i.bak 's/.fas//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
			fi
		else
			#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
			sed -i.bak 's/(/_/g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
			sed -i.bak 's/)//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
			#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
			sed -i.bak 's/_contigs//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
			sed -i.bak 's/.fas//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
		fi
	done
	
	#Prepare concatenated dataset and transform it to phylip format
	if [[ $AMAS =~ "slow" ]]; then
		#Much slower option but works also in case of many genes
		xx=0
		for f in *.fas ; do
			echo $f
			if [ $xx -eq 0 ]; then
				cp $f concatenated.workfasta
				xx=$((xx + 1))
			else
				python3 AMAS.py concat -i concatenated.workfasta $f -f fasta -d dna -u fasta -t concatenated.workfasta2 >/dev/null
				mv concatenated.workfasta2 concatenated.workfasta
				xx=$((xx + 1))
			fi
		done
		mv concatenated.workfasta concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta
		python3 AMAS.py convert -d dna -f fasta -i concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta -u phylip >/dev/null
		mv concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta-out.phy concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip
	else
		#Faster solution but with really many genes generate 'Argument list too long' error
		python3 AMAS.py concat -i *.fas -f fasta -d dna -u fasta -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta >/dev/null
		python3 AMAS.py concat -i *.fas -f fasta -d dna -u phylip -t concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip >/dev/null
	fi
	#Calculate proportion of missing data in the concatenated alignment
	#Removes line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > tmp && mv tmp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta
	#Calculate length of alignment: 1. get second line and count length, 2. decrease value by one (because previous command also counted LF)
	length=$(cat concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta | head -n 2 | tail -n 1 | wc -c)
	length=`expr $length - 1`
	#Replace newline with ' ' if line starts with '>' (i.e., merge headers with data into single line separated by space)
	cat concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta | sed '/^>/{N; s/\n/ /;}' > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.modif.fasta
	#Cut first part until space, i.e. header, and remove '>'
	cat concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.modif.fasta | cut -f1 -d" " | sed 's/>//' > headers.txt
	#Cut only part after the first space, i.e., only sequence, change all missing data (-, ?, N) to 'n', replace all other characters then 'n' by nothing and print percentage of 'n's in each sequence
	cat concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.modif.fasta | cut -f2 -d" " | sed 's/[?N-]/n/g' | sed 's/[^n]//g' | awk -v val=$length '{ print (length*100)/val }' > missingpercentage.txt
	paste headers.txt missingpercentage.txt > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt
	#Calculate mean of all values
	echo -e "MEAN\t$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt)" > mean.txt
	cat concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt mean.txt > tmp && mv tmp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt
	rm headers.txt missingpercentage.txt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.modif.fasta mean.txt
	
	#Modify partition file: take 2nd and following parts of each line separated by "_" (remove 'p1_' etc. introduced by AMAS) | add ';' at the end of each line
	cat partitions.txt | cut -d"_" -f2- | awk '{ print $0 ";" }' | sed 's/-/ - /g' | sed 's/=/ = /g' > part.file
	#Copy concatenated file to home
	if [[ $update =~ "yes" ]]; then
		if [[ $requisite =~ "yes" ]]; then
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.fasta
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.phylip
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/
		else
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
		fi
	else
		if [[ $requisite =~ "yes" ]]; then
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.fasta
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_with_requisite.phylip
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
		else
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
			cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
		fi
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
	#wget --no-check-certificate https://github.com/brettc/partitionfinder/archive/v2.1.1.tar.gz
	#tar xfz v2.1.1.tar.gz
	cp $source/partitionfinder-2.1.1.tar.gz .
	tar xfz partitionfinder-2.1.1.tar.gz
	mv partitionfinder-2.1.1 partitionfinder
	#cp $source/partitionfinder-2.0.0-pre13.tar.gz .
	#tar xfz partitionfinder-2.0.0-pre13.tar.gz
	#mv partitionfinder-2.0.0-pre13 partitionfinder
	
	if [[ $PBS_O_HOST == *".cz" ]]; then
		#Remove python3 module and dd python modules necessary for PartitionFinder
		module rm python-3.4.1-gcc
		#module add python-2.7.10-gcc
		module add python27-modules-gcc
		module add hdf5-1.8.12-gcc
	fi
	
	# Run PartitionFinder
	#python partitionfinder/PartitionFinder.py PartitionFinder --raxml --rcluster-max 1000 --rcluster-percent 10
	if [[ $PBS_O_HOST == *".cz" ]]; then
		python partitionfinder/PartitionFinder.py -p $TORQUE_RESC_TOTAL_PROCS PartitionFinder --raxml --no-ml-tree
	else
		python partitionfinder/PartitionFinder.py -p $numbcores PartitionFinder --raxml --no-ml-tree
	fi
	#Prepare RAxML (ExaML) partition file from PartitionFinder results
	grep "DNA, Subset" PartitionFinder/analysis/best_scheme.txt > RAxMLpartitions.txt
	echo -e "\nPartitionFinder finished...\n"
	#Copy results to home
	if [[ $update =~ "yes" ]]; then
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
		cp -r PartitionFinder $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
		cp RAxMLpartitions.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	else
		mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
		cp -r PartitionFinder $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
		cp RAxMLpartitions.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	fi
fi

#Check if there are completely undetermined columns in alignment (RAxML -y will produced .reduced alignment and partition files)
#Compute parsimony tree only and produce reduced alignment and appropriate reduced partition file
echo -e "Computing parsimony tree to check undetermined columns...\n"
$raxmlseq -y -m GTRCAT -p 12345 -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -q RAxMLpartitions.txt -n check > RAxML_to_check.log
#Test if reduced files were produced, copy them to home and rename to original names
if [ -f RAxMLpartitions.txt.reduced ]; then
	echo "Reduced alignment found...using it"
	if [[ $update =~ "yes" ]]; then
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.reduced $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
		cp RAxMLpartitions.txt.reduced $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	else
		cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.reduced $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
		cp RAxMLpartitions.txt.reduced $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
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
parse-examl -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -m DNA -q RAxMLpartitions.txt -n concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.partitioned >> parse-examl_besttree.log
#Prepare starting tree using RAxML
echo -e "Preparing starting tree using RAxML...\n"
$raxmlseq -y -m GTRCAT -p 12345 -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -n RandomStartingTree > RAxML_starttree_besttree.log
#Run ExaML
echo -e "Building best ExaML tree...\n"
if [[ $PBS_O_HOST == *".cz" ]]; then
	mpirun $examlbin -t RAxML_parsimonyTree.RandomStartingTree -m GAMMA -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.partitioned.binary -n examl.tre > ExaMLbest.log 2>>ExaML_error.log
else
	$examlbin -t RAxML_parsimonyTree.RandomStartingTree -m GAMMA -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.partitioned.binary -n examl.tre > ExaMLbest.log 2>>ExaML_error.log
fi
#Measure time
timenow=`date +%s`
timespan=`echo "scale=2;($timenow-$time)/60" | bc`

#Delete checkpoint files
rm *Check* 2>/dev/null

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp ExaML* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	echo -e "BestML done...in $timespan mins" >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaML_bootstrap_progress.txt
else
	cp ExaML* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	echo -e "BestML done...in $timespan mins" >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaML_bootstrap_progress.txt
fi

if [[ $examlboot =~ "yes" ]]; then
	#Bootstrap ExaML
	#generate 100 bootstrap replicates
	echo -e "Preparing bootstrap replicates using RAxML...\n"
	$raxmlseq -N 100 -b 12345 -f j -m GTRCAT -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip -q RAxMLpartitions.txt -n REPS > RAxML_preparingBS.log
	
	for (( i=0; i<=99; i++ ))
	do
		echo -en "Bootstrap replicate $i..."\\r
		#Set starting time
		time=`date +%s`
		#generate parsimony starting tree
		$raxmlseq -y -m GTRCAT -p 12345 -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i} -q RAxMLpartitions.txt.BS${i} -n T${i} >> RAxML_boot.log
		#Prepare binary file for ExaML
		parse-examl -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i} -m DNA -q RAxMLpartitions.txt.BS${i} -n concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i} >> parse-examl_boot.log
		#Run ExaML
		if [[ $PBS_O_HOST == *".cz" ]]; then
			mpirun $examlbin -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i}.binary -m GAMMA -t RAxML_parsimonyTree.T${i} -n BINF_${i} >> ExaML_boot.log 2>>ExaML_error.log
		else
			$examlbin -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i}.binary -m GAMMA -t RAxML_parsimonyTree.T${i} -n BINF_${i} >> ExaML_boot.log 2>>ExaML_error.log
		fi
		mpirun $examlbin -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip.BS${i}.binary -m GAMMA -t RAxML_parsimonyTree.T${i} -n BINF_${i} >> ExaML_boot.log 2>>ExaML_error.log
		#Measure time
		timenow=`date +%s`
		timespan=`echo "scale=2;($timenow-$time)/60" | bc`
		if [[ $update =~ "yes" ]]; then
			echo -e "Bootstrap $i done...in $timespan mins" >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaML_bootstrap_progress.txt
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaMLfiles
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaMLfiles/bootstrap_trees
			cp ExaML_result.BINF_${i} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaMLfiles/bootstrap_trees
			
		else
			echo -e "Bootstrap $i done...in $timespan mins" >> $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaML_bootstrap_progress.txt
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaMLfiles
			mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaMLfiles/bootstrap_trees
			cp ExaML_result.BINF_${i} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaMLfiles/bootstrap_trees
		fi
	done
	
	#Combine all bootstrap trees
	echo -e "Combining bootstrap trees...\n"
	cat ExaML_result.BINF* > ExaML_bootstrap.tre
	#Map bootstrap values onto the single best tree
	echo -e "Mapping bootstrap support values onto bestML tree...\n"
	$raxmlseq -f b -m GTRGAMMA -z ExaML_bootstrap.tre -t ExaML_result.examl.tre -n ExaML_Bootstrap > RAxML_mapBStoBestML.log
fi

#Delete checkpoint files
rm *Check* 2>/dev/null
#Delete bootstrap trees
rm ExaML_result.BINF* 2>/dev/null

#Rename some files
mv ExaML_result.examl.tre ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre
if [[ $examlboot =~ "yes" ]]; then
	mv ExaML_bootstrap.tre ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.trees
	mv RAxML_bipartitions.ExaML_Bootstrap ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre
fi

#(Re)root a final concatenated species trees with $OUTGROUP
if [ -n "$OUTGROUP" ]; then
	echo -e "Rooting trees...\n"
	nw_reroot -s ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre $OUTGROUP > tmp && mv tmp ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre
	nw_reroot -s ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.trees $OUTGROUP > tmp && mv tmp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.trees
	nw_reroot -s ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre $OUTGROUP > tmp && mv tmp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre
fi

#Copy results to home
if [[ $update =~ "yes" ]]; then
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaMLfiles
	cp *ExaML* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML/ExaMLfiles
	cp ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	if [[ $examlboot =~ "yes" ]]; then
		cp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.trees $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
		cp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
	fi
else
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaMLfiles
	cp *ExaML* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML/ExaMLfiles
	cp ExaML_BestML_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	if [[ $examlboot =~ "yes" ]]; then
		cp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.trees $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
		cp ExaML_bootstrap_${MISSINGPERCENT}_${SPECIESPRESENCE}${modif1}.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
	fi
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/${modif}concatenatedExaML
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/${modif}concatenatedExaML
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08f
fi

echo -e "\nScript HybPhyloMaker 8f finished...\n"
