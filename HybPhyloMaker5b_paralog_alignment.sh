#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker5b_paralog_alignment_and_FastTree_for_selected
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker5b_paralog_alignment_and_FastTree_for_selected
#$ -o HybPhyloMaker5b_paralog_alignment_and_FastTree_for_selected.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                 Script 05b - MAFFT and FastTree for paralogs                 *
# *                                   v.1.8.0b                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Make alignments of loci and theirs paralogs and makes a common OP trees using FastTree
# Selection is based on maximum missing data per sample allowed ($MISSINGPERCENT) and minimum species percentage presence per assembly ($SPECIESPRESENCE)
# Edit these values in settings.cfg
# Run first HybPhyloMaker5_missingdataremoval.sh with the same settings
# Only works for datasets treated with Paralog Wizard !!!

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker5b is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add fasttree-2.1.8
	module add mafft-7.029
	module add debian8-compat #for proper working of 'python3' on some computing nodes
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker5b is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir05b
	cd workdir05b
	#Add necessary modules
	module load bioinformatics/fasttree/2.1.8
	module load bioinformatics/raxml/8.2.11
	module load tools/R/3.4.1
	module load bioinformatics/anaconda3/5.1 #for python3
else
	echo -e "\nHybPhyloMaker5b is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir05b
	cd workdir05b
fi

#Write log
logname=HPM5b
echo -e "HybPhyloMaker5b: MAFFT and FastTree for paralogs after ParalogWizard" > ${logname}.log
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
for set in data cp corrected MISSINGPERCENT SPECIESPRESENCE FastTreeBoot; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
else
	echo -e "Working with exons"
	type="exons"
fi

#Settings for (un)corrected reading frame
if [[ $corrected =~ "yes" ]]; then
	alnpath=$type/80concatenated_exon_alignments_corrected
	alnpathselected=$type/81selected_corrected
	treepath=$type/82trees_corrected
	echo -e "...with corrected reading frame\n"
else
	alnpath=$type/70concatenated_exon_alignments
	alnpathselected=$type/71selected
	treepath=$type/72trees
	echo -e "\n"
fi

#Check necessary file
echo -ne "Testing if input data are available..."
if [ -f "$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt" ]; then
	if [ -d "$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}" ]; then
		if [ "$(ls -A $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT})" ]; then
			echo -e "OK\n"
		else
			echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is empty. Exiting...\n"
			rm -d ../workdir05b/ 2>/dev/null
			exit 3
		fi
	else
		echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}' is missing. Exiting...\n"
		rm -d ../workdir05b/ 2>/dev/null
		exit 3
	fi
else
	echo -e "'$path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt' is missing. Exiting...\n"
	rm -d ../workdir05b/ 2>/dev/null
	exit 3
fi

#Test if folder for results exits
if [ -d "$path/${alnpathselected}${MISSINGPERCENT}/OPtrees" ]; then
	echo -e "Directory '$path/${treepath}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees' already exists. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir05b/ 2>/dev/null
	exit 3
else
	if [[ ! $location == "1" ]]; then
		if [ "$(ls -A ../workdir05b)" ]; then
			echo -e "Directory 'workdir05b' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
			rm -d ../workdir05b/ 2>/dev/null
			exit 3
		fi
	fi
fi


#Add necessary scripts and files
cp $path/${alnpathselected}${MISSINGPERCENT}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt .
#Prepare list of paralogous loci
grep "para" selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d '_' -f2 | sed 's/para//' > list.txt
grep -f list.txt selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt > tmp && mv tmp selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
# Copy and modify selected FASTA files
echo -en "\nCopying and modifying selected FASTA files..."
for i in $(cat selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt)
do
	cp $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas .
	#Substitute '(' by '_' and ')' by nothing ('(' and ')' not allowed in FastTree)
	sed -i.bak 's/(/_/g' *${i}_modif${MISSINGPERCENT}.fas
	sed -i.bak 's/)//g' *${i}_modif${MISSINGPERCENT}.fas
	#Delete '_contigs' and '.fas' from labels (i.e., keep only genus-species_nr)
	sed -i.bak 's/_contigs//g' *${i}_modif${MISSINGPERCENT}.fas
	sed -i.bak 's/.fas//g' *${i}_modif${MISSINGPERCENT}.fas
done
#Make a list of all fasta files
ls *.fas | cut -d"." -f1 > FileForOPTree.txt
#Make dir for results
mkdir -p $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees

echo -e "done\n"

#----------------Generate alignments of locus and its paralogue----------------
echo -e "Name\tLocus\tParalocus" > parastat.txt
cp parastat.txt $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees

for file in $(cat list.txt); do
	#Calculate number of samples
	locus=$(grep ">" *${file}_modif70.fas | wc -l)
	paralocus=$(grep ">" *${file}para_modif70.fas | wc -l)
	echo -e "${file}\t${locus}\t${paralocus}" >> parastat.txt
	cp parastat.txt $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees
	#Rename paralogs
	sed -i 's/>/>para_/' *${file}para_modif70.fas
	#Merge locus and paralocus alignments
	cat *${file}_modif70.fas *${file}para_modif70.fas > ${file}ANDpara_modif70.fas
	#Alignment using MAFFT
	mafft --auto ${file}ANDpara_modif70.fas > ${file}ANDpara_modif70.mafft
	cp ${file}ANDpara_modif70.mafft $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees
	#Tree with FastTree
	fasttreemp -nt ${file}ANDpara_modif70.mafft > ${file}ANDpara_modif70.fast.tre
	cp ${file}ANDpara_modif70.fast.tre $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees
done

#Combine all trees
cat *fast.tre > allOPtrees.tre
cp allOPtrees.tre $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/${alnpathselected}${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}/OPtrees

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir05b
fi

exit

#----------------Generate gene trees using FastTree----------------
echo -e "Generating FastTrees..."
numbertrees=$(cat FileForFastTree.txt | wc -l)
calculating=0
for file in $(cat FileForFastTree.txt); do
	calculating=$((calculating + 1))
	#FastTree
	echo -e "$file ($calculating out of $numbertrees)"
	echo -e "\nProcessing file: ${file}" >> FastTree.log
	if [[ $location == "1" ]]; then
		fasttreemp -nt ${file}.fas > ${file}.fast.tre 2>>FastTree.log
	elif [[ $location == "2" ]]; then
		FastTreeMP -nt ${file}.fas > ${file}.fast.tre 2>>FastTree.log
	else
		$fasttreebin -nt ${file}.fas > ${file}.fast.tre 2>>FastTree.log
	fi
	cp *$file.fast.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
done
#Copy log to home
cp FastTree.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree

#----------------Generate bootstrapped gene trees using FastTree----------------
#Bootstrap using FastTree 
if [[ $FastTreeBoot =~ "yes" ]]; then
	echo -e "\nGenerating bootstrapped FastTrees..."
	numbertrees=$(cat FileForFastTree.txt | wc -l)
	calculating=0
	for file in $(cat FileForFastTree.txt); do
		calculating=$((calculating + 1))
		#Generate 100 replicated datasets using RAxML -f j
		echo -e "\n$file ($calculating out of $numbertrees)"
		echo -e "\nProcessing file: ${file}" >> FastTreeBoot.log
		echo -e "\nGenerating bootstrap replicates for: ${file}" >> raxml.log
		if [[ $location == "1" ]]; then
			raxmlHPC -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS >> raxml.log
		elif [[ $location == "2" ]]; then
			raxmlHPC-SSE3 -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS >> raxml.log
		else
			echo -ne "Generating BS replicates..."
			$raxmlseq -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS >> raxml.log
		fi
		#Test if there was an error in generating BS replicate
		if grep -Fq "Error in BS replicate" RAxML_info.BS; then
			echo -e "warning: reduced alignment used"
			python3 AMAS.py convert -d dna -f phylip -i ${file}.fas.reduced -u fasta >/dev/null
			mv ${file}.fas.reduced-out.fas ${file}.fas
			rm RAxML_info.BS
			echo -e "\nGenerating bootstrap replicates for: ${file}.reduced" >> raxml.log
			if [[ $location == "1" ]]; then
				raxmlHPC -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS >> raxml.log
			elif [[ $location == "2" ]]; then
				raxmlHPC-SSE3 -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS >> raxml.log
			else
				$raxmlseq -f j -b 12345 -N 100 -s ${file}.fas -m GTRCAT -n BS >> raxml.log
			fi
		else
			echo "OK"
		fi
		python3 AMAS.py convert -d dna -f phylip -i *fas.BS* -u fasta >/dev/null
		#Loop over replicates and calculate FastTree for each of them
		for i in {0..99}; do
			echo -en "Replicate $i"\\r
			echo -e "\nFile: $file, replicate $i" >>FastTreeBoot.log
			if [[ $location == "1" ]]; then
				fasttreemp -nt ${file}.fas.BS${i}-out.fas > ${file}.BS${i}.fast.tre 2>>FastTreeBoot.log
			elif [[ $location == "2" ]]; then
				FastTreeMP -nt ${file}.fas.BS${i}-out.fas > ${file}.BS${i}.fast.tre 2>>FastTreeBoot.log
			else
				$fasttreebin -nt ${file}.fas.BS${i}-out.fas > ${file}.BS${i}.fast.tre 2>>FastTreeBoot.log
			fi
		done
		#Combine all bootstrap trees to a single file
		cat ${file}.BS*.fast.tre > ${file}.boot.fast.trees
		#Delete RAxML_info.BS, BS files and trees and BS alignments in fasta format (generated by AMAS)
		rm ${file}.fas.BS*
		rm -f RAxML_info.BS
		rm ${file}.BS*.fast.tre
		#Map bootstrap support values onto the original tree
		perl -I . ./CompareToBootstrap.pl -tree *${file}.fast.tre -boot ${file}.boot.fast.trees > ${file}.boot.fast.tre
		#Copy bootstrap trees and final tree with bootstrap values to home
		cp ${file}.boot.fast.trees $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
		cp ${file}.boot.fast.tre $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	done
	#Copy logs to home
	cp *.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
fi

#----------------Make a summary table with statistical properties for trees using R----------------
#Remove 'debian8-compat' module if on MetaCentrum, otherwise R packages are not correctly loaded
if [[ $location == "1" ]]; then
	module rm debian8-compat
	module add r/4.4.0-gcc-10.2.1-ssuwpvb
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"
fi
#Copy script
cp $source/tree_props.r .
cp $source/treepropsPlot.r .
cp $source/LBscores.R .
mkdir trees
mkdir alignments
#Copy all fasta alignments to subfolder 'alignments'
cp *.fas alignments/
#Copy all tree files (*.tre) but not bootstrap trees (*boot*) to subfolder 'trees'
cp `find ./*.tre -maxdepth 1 ! -name '*boot*'` trees/
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
cp LBscores.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
cp LBscoresPerTaxon.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
cp LBscoresSDPerLocus.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
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
cp tree_stats_table.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
cp *.png $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
#Remove results
rm *.png
rm LB*.txt LB*.csv
mv tree_stats_table.csv tree_stats_table_noboot.csv

#----------------Make a summary table for bootstrapped trees----------------
if [[ $FastTreeBoot =~ "yes" ]]; then
	#Remove everything from subfolder 'trees'
	rm ./trees/*.*
	#Copy all trees with bootstrap values to subfolder 'trees'
	cp `find ./*boot.fast.tre -maxdepth 1` trees/
	#Run R script for tree properties calculation
	echo -e "Calculating tree properties for bootstrapped trees...\n"
	echo -e "\nCalculating tree properties for bootstrapped trees using tree_props.r\n" >> R.log
	R --slave -f tree_props.r >> R.log 2>&1
	#Run R script for calculation of LB score
	echo -e "Calculating and parsing LB score...\n"
	echo -e "\nCalculating and parsing LB score for bootstrapped trees using LBscores.R\n" >> R.log
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
	mv LBscores.csv LBscores_bootstrap.csv
	mv LBscoresPerTaxon.txt LBscoresPerTaxon_bootstrap.txt
	mv LBscoresSDPerLocus.txt LBscoresSDPerLocus_bootstrap.txt
	cp LBscores_bootstrap.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	cp LBscoresPerTaxon_bootstrap.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	cp LBscoresSDPerLocus_bootstrap.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	#Combine 'LBscoresSDPerLocus_bootstrap.txt' with 'tree_stats_table.csv'
	awk '{ print $2 }' LBscoresSDPerLocus_bootstrap.txt > tmp && mv tmp LBscoresSDPerLocus_bootstrap.txt
	paste tree_stats_table.csv LBscoresSDPerLocus_bootstrap.txt | tr "\t" "," > tmp && mv tmp tree_stats_table.csv
	#Replace 'NaN' by '0' (otherwise following plotting in R will not work)
	sed -i.bak 's/NaN/0/g' tree_stats_table.csv
	echo -e "Plotting boxplots/histograms for bootstrapped tree properties...\n"
	echo -e "\nPlotting boxplots/histograms for bootstrapped tree properties using treepropsPlot.r\n" >> R.log
	if [[ $location == "1" ]]; then
		#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
		#xvfb-run R --slave -f treepropsPlot.r
		R --slave -f treepropsPlot.r >> R.log 2>&1
	else
		R --slave -f treepropsPlot.r >> R.log 2>&1
	fi
	#Rename results and copy results to home
	mv tree_stats_table.csv tree_stats_table_bootstrap.csv
	cp tree_stats_table_bootstrap.csv $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	#Rename all PNG files generated by R (add 'bootstrap')and copy them to home
	for file in *.png; do mv "$file" "${file/.png/_bootstrap.png}"; done
	cp *.png $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
fi

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
sed -i.bak 's/,/\t/g' tree_stats_table_noboot.csv
#If FasTree bootstrap was not requested, create "fictive" empty tree_stats_table_bootstrap.csv (to avoid error messages on next lines)
#This "fictive" empty tree_stats_table_bootstrap.csv will be deleted 
if [[ $FastTreeBoot =~ "no" ]]; then
	touch tree_stats_table_bootstrap.csv
fi
sed -i.bak 's/,/\t/g' tree_stats_table_bootstrap.csv
#Take first line as a header
head -n1 tree_stats_table_noboot.csv > head2.txt
#Remove first line | sort
sed 1d tree_stats_table_noboot.csv | sort > tree_stats_table_noboot_sorted.csv
sed 1d tree_stats_table_bootstrap.csv | sort > tree_stats_table_bootstrap_sorted.csv
#Combine header and sorted file
cat head2.txt tree_stats_table_noboot_sorted.csv > tmp && mv tmp tree_stats_table_noboot_sorted.csv
cat head2.txt tree_stats_table_bootstrap_sorted.csv > tmp && mv tmp tree_stats_table_bootstrap_sorted.csv
#Combine both files
paste summarySELECTED_sorted.txt tree_stats_table_noboot_sorted.csv > combined_noboot.txt
paste summarySELECTED_sorted.txt tree_stats_table_bootstrap_sorted.csv > combined_bootstrap.txt
#Delete bootstrap-related files if no bootstrap requested
if [[ $FastTreeBoot =~ "no" ]]; then
	rm *bootstrap*
fi
#Rename colums
sed -i.bak 's/Alignment_length/Aln_length/' combined_*.txt
sed -i.bak 's/Missing_percent/Missing_perc/' combined_*.txt
sed -i.bak 's/Proportion_parsimony_informative/Prop_pars_inf/' combined_*.txt
sed -i.bak 's/MstatX_entropy/Aln_entropy/' combined_*.txt
sed -i.bak 's/Average_bootstrap/Bootstrap/' combined_*.txt
sed -i.bak 's/Average_branch_length/Branch_length/' combined_*.txt
sed -i.bak 's/Avg_p_dist/P_distance/' combined_*.txt
sed -i.bak 's/Slope/Satur_slope/' combined_*.txt
sed -i.bak 's/R_squared/Satur_R_sq/' combined_*.txt
sed -i.bak 's/LBscoreSD/LBscore_SD/' combined_*.txt

#Run comparison plots for FastTrees without true bootstrap ('noboot')
mv combined_noboot.txt combined.txt
echo -e "Plotting gene properties correlations for FastTrees without true bootstrap...\n"
echo -e "\nPlotting gene properties correlations for FastTrees without true bootstrap using plotting_correlations.R\n" >> R.log
if [[ $location == "1" ]]; then
	#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
	#xvfb-run R --slave -f plotting_correlations.R
	R --slave -f plotting_correlations.R >> R.log 2>&1
else
	R --slave -f plotting_correlations.R >> R.log 2>&1
fi
cp genes_corrs.* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
rm genes_corrs.*
mv combined.txt gene_properties.txt
cp gene_properties.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree

#Run comparison plots for FastTrees with true bootstrap ('bootstrap')
if [[ $FastTreeBoot =~ "yes" ]]; then
	echo -e "Plotting gene properties correlations for bootstrapped FastTrees...\n"
	echo -e "\nPlotting gene properties correlations for bootstrapped FastTrees using plotting_correlations.R\n" >> R.log
	mv combined_bootstrap.txt combined.txt
	if [[ $location == "1" ]]; then
		#Run R script for boxplot/histogram visualization (run via xvfb-run to enable generating PNG files without X11 server)
		#xvfb-run R --slave -f plotting_correlations.R
		R --slave -f plotting_correlations.R >> R.log 2>&1
	else
		R --slave -f plotting_correlations.R >> R.log 2>&1
	fi
	mv genes_corrs.png genes_corrs_bootstrap.png
	mv genes_corrs.pdf genes_corrs_bootstrap.pdf
	cp genes_corrs_bootstrap.* $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
	rm genes_corrs_bootstrap.*
	mv combined.txt gene_properties_bootstrap.txt
	cp gene_properties_bootstrap.txt $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree
fi

#Copy R.log to home
cp R.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/FastTree

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir06b
fi

echo -e "HybPhyloMaker 5b finished...\n"
