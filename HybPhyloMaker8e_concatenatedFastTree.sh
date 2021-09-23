#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=24gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8e_concatenated_tree
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8e_concatenated_tree
#$ -o HybPhyloMaker8e_concatenated_tree.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *            Script 08e - concatenated species tree using FastTree             *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Compute species tree from selected concatenated genes
#Take genes specified in /71selected${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
#from /71selected/deleted_above${MISSINGPERCENT}
#Run first
#(1) HybPhyloMaker5_missingdataremoval.sh with the same ${MISSINGPERCENT} and ${SPECIESPRESENCE} values

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8e is running on MetaCentrum..."
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
	module add newick-utils-13042016
	module add debian8-compat
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8e is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08e
	cd workdir08e
	#Add necessary modules
	module load bioinformatics/fasttree/2.1.8
	module load bioinformatics/anaconda3/5.1 #python3 and NewickUtilities
else
	echo -e "\nHybPhyloMaker8e is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08e
	cd workdir08e
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

#Check necessary file
echo -ne "\nTesting if input data are available..."
if [[ $update =~ "yes" ]]; then
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt" ]; then
		#if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
		if [ 0 -lt $(find $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT} -maxdepth 1 -name "*ssembly_*_modif${MISSINGPERCENT}.fas" -exec ls {} + 2>/dev/null | wc -w) ]; then #to avoid 'Argument list too long' error
			echo -e "OK\n"
		else
			echo -e "no alignmenet files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
			rm -d ../workdir08e 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/updatedSelectedGenes/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}_update.txt' is missing. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
else
	if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
		#if [ 0 -lt $(ls $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_*_modif${MISSINGPERCENT}.fas 2>/dev/null | wc -w) ]; then
		if [ 0 -lt $(find $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT} -maxdepth 1 -name "*ssembly_*_modif${MISSINGPERCENT}.fas" -exec ls {} + 2>/dev/null | wc -w) ]; then #to avoid 'Argument list too long' error
			echo -e "OK\n"
		else
			echo -e "no alignmenet files in FASTA format found in '$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}'. Exiting..."
			rm -d ../workdir08e 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
fi

#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
else
	if [ -d "$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08e)" ]; then
		echo -e "Directory 'workdir08e' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08e 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8e
echo -e "HybPhyloMaker8e: concatenated species tree using FastTree" > ${logname}.log
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
for set in data selection cp corrected update MISSINGPERCENT SPECIESPRESENCE tree OUTGROUP; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

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
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	mkdir -p $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

# Copy and modify selected FASTA files
echo -e "Copying and modifying FASTA files...\n"
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d"_" -f2); do
	#If working with 'corrected' copy trees starting with 'CorrectedAssembly'
	#cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/*ssembly_${i}_modif${MISSINGPERCENT}.fas .
	find $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT} -maxdepth 1 -name "*ssembly_${i}_modif${MISSINGPERCENT}.fas" -exec cp -t . {} + #to avoid 'Argument list too long' error
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in RAxML/FastTree)
	#sed -i.bak 's/(/_/g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	find . -maxdepth 1 -name "*ssembly_${i}_modif${MISSINGPERCENT}.fas" | xargs sed -i.bak 's/(/_/g' #to avoid 'Argument list too long' error
	#sed -i.bak 's/)//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	find . -maxdepth 1 -name "*ssembly_${i}_modif${MISSINGPERCENT}.fas" | xargs sed -i.bak 's/)//g' #to avoid 'Argument list too long' error
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	#sed -i.bak 's/_contigs//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	find . -maxdepth 1 -name "*ssembly_${i}_modif${MISSINGPERCENT}.fas" | xargs sed -i.bak 's/_contigs//g' #to avoid 'Argument list too long' error
	#sed -i.bak 's/.fas//g' *ssembly_${i}_modif${MISSINGPERCENT}.fas
	find . -maxdepth 1 -name "*ssembly_${i}_modif${MISSINGPERCENT}.fas" | xargs sed -i.bak 's/.fas//g' #to avoid 'Argument list too long' error
done

#Prepare concatenated dataset and transform it to phylip format
echo -e "Concatenating...\n"
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
echo -e "Calculating proportion of missing data...\n"
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

#Copy concatenated file to home
if [[ $update =~ "yes" ]]; then
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.phylip $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}_missingperc.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

#FastTree for concatenated dataset
echo -e "Computing FastTree for concatenated dataset...\n"
if [[ $location == "1" ]]; then
	fasttreemp -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
elif [[ $location == "2" ]]; then
	FastTreeMP -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
else
	$fasttreebin -nt concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fasta > concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
fi

#Removing '_cpDNA' from names in alignment
sed -i.bak 's/_cpDNA//g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre

#(Re)root a final concatenated species tree with $OUTGROUP
if [ -n "$OUTGROUP" ]; then
	nw_reroot -s concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $OUTGROUP > tmp && mv tmp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
fi

#Modify labels in concatenated trees
sed -i.bak2 's/-/ /g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre
sed -i.bak2 's/_/ /g' concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	cp concatenated${MISSINGPERCENT}_${SPECIESPRESENCE}.fast.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/concatenated
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/concatenated
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08e
fi

echo -e "\nHybPhyloMaker8e finished...\n"
