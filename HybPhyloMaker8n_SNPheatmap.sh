#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8n_SNPheatmap
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8n_SNPheatmap
#$ -o HybPhyloMaker8n_SNPheatmap.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                           Script 08n - SNP heatmap                           *
# *                                   v.1.8.0b                                   *
# *                                  Tomas Fer                                   *
# *       Dept. of Botany, Charles University, Prague, Czech Republic, 2025      *
# *                           tomas.fer@natur.cuni.cz                            *
# ********************************************************************************

#Plot heatmap based on all filtered SNPs from concatenated alignment
#Take ALL gene alignments, concatenate them using AMAS and creates VCF file of variable sites using SNP-sites

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8n is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add newick-utils-13042016
	module add mambaforge #to add SNP-sites later
	module add bcftools-1.11
	module add vcftools-0.1.16
	module unload python-2.6.6-gcc #avoid conflict with the next module
	#module add python36-modules-gcc
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8n is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08n
	cd workdir08n
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #adds NewickUtilities
	module load ??? #add SNP-sites on Hydra (maybe via conda?)
	module add bcftools???
	module add vcftools???
else
	echo -e "\nHybPhyloMaker8n is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08n
	cd workdir08n
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
echo -e "\nTesting if input data are available..."
#test for folder with alignments
#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNPheatmap/${snpmiss}missing" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNPheatmap/${snpmiss}missing' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08n 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNPheatmap/${snpmiss}missing" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNPheatmap/${snpmiss}missing' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08n 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08n)" ]; then
		echo -e "Directory 'workdir08n' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08n 2>/dev/null
		exit 3
	fi
fi

#Check if 'snpmiss' is integer
if [[ ! "$snpmiss" =~ ^[0-9]+$ ]]; then
	echo -e "Variable 'snpmiss' is not correctly set:\nNow: $snpmiss\nShould be an integer between 0 and 100.\n Exiting...\n"
	rm -d ../workdir08n 2>/dev/null
	exit 3
else
	snpmiss=$((10#${snpmiss}))
	if [[ "$snpmiss" -gt 100 ]]; then
		echo -e "Variable 'snpmiss' is not correctly set:\nNow: $snpmiss\nShould be an integer between 0 and 100.\n Exiting...\n"
		rm -d ../workdir08n 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8n
echo -e "HybPhyloMaker8n: SNP heatmap" > ${logname}.log
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
for set in data selection cp corrected update tree OUTGROUP SNPs snpmiss; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done
if [ ! -z "$selection" ]; then
	echo -e "\nList of excluded samples" >> ${logname}.log
	cat $source/excludelist.txt >> ${logname}.log
	echo >> ${logname}.log
fi

#Add necessary scripts and files
echo -e "\nCopying scripts...\n"
cp $source/AMAS.py .
cp $source/SNPheatmap.R .

#Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNPheatmap
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNPheatmap
fi

#Concatenate gene alignments and create VCF
echo -e "Copying and modifying FASTA files...\n"
#copy all FASTA files (not selected only!)
find $path/${alnpath}/ -maxdepth 1 -name "*ssembly_*.fasta" -exec cp -t . {} + #to avoid 'Argument list too long' error
#concatenate FASTA files
echo -e "Concatenating...\n"
if [[ $AMAS =~ "slow" ]]; then
	#Much slower option but works also in case of many genes
	xx=0
	for f in *.fasta ; do
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
	mv concatenated.workfasta concatenated.fasta
else
	#Faster solution but with really many genes generate 'Argument list too long' error
	python3 AMAS.py concat -i *.fasta -f fasta -d dna -u fasta -t concatenated.fasta >/dev/null
fi
rm *Assembly*.fasta

echo -e "Modifying concatenated FASTA...\n"
#remove line breaks from FASTA
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' concatenated.fasta > tmp && mv tmp concatenated.fasta
sed -i 's/_contigs.fas//' concatenated.fasta
#replace '?'s by 'N'
sed -i 's/\?/N/g' concatenated.fasta
#replace 'N's by '-'
sed -i '/^>/!s/N/-/g' concatenated.fasta

if [[ $snpmiss -eq 0 ]]; then
	echo -e "Running standard SNP-site...\n"
	if [[ $PBS_O_HOST == *".cz" ]]; then
		mamba activate SNP-sites
		snp-sites -v concatenated.fasta > concatenated.vcf
		mamba deactivate
		module unload mambaforge
	else
		snp-sites -v concatenated.fasta > concatenated.vcf
	fi
	#Filter VCF
	echo -e "Filtering VCF...\n"
	gzip concatenated.vcf
	bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2' -m2 -M2 -O z -o concatenated_bcf.vcf.gz concatenated.vcf.gz
	gzip -d -c concatenated_bcf.vcf.gz > concatenated_bcf.vcf
	mv concatenated_bcf.vcf.gz FILTEREDconcatenated_${snpmiss}missing.vcf.gz
	mv concatenated.vcf.gz FULLconcatenated_${snpmiss}missing.vcf.gz
	mv concatenated_bcf.vcf FILTEREDconcatenated_${snpmiss}missing.vcf
	mkdir 
else
	#change alignment gaps ('-') to 'X'
	echo -e "Running SNP-site with gaps changed to 'X'...\n"
	sed '/^>/! s/-/X/g' concatenated.fasta > concatenatedModif.fasta
	if [[ $PBS_O_HOST == *".cz" ]]; then
		mamba activate SNP-sites
		snp-sites -v concatenatedModif.fasta > concatenatedModif.vcf
		mamba deactivate
		module unload mambaforge
	else
		snp-sites -v concatenatedModif.fasta > concatenatedModif.vcf
	fi
	#separate header and data
	grep "^#" concatenatedModif.vcf > header.txt
	grep -v "^#" concatenatedModif.vcf > data.txt
	#remove triallelic/tetraallelic SNPs
	grep -v "[ACTGX],[ACTGX],[ACTGX]" data.txt > dataToModif.txt
	#Modification commands - #explanation:
	#-F"\" - works on TABs as column separators
	#always checks value in 5th colum ($5), i.e. ALT
	#if $5=X then set the value to "*" (i.e. indel) and then checks the values in columns from 10 until the end (for (i=10; i<=NF; i++)) - i.e. genotype values - and change all "1"s to "." (i.e. missing)
	#than it checks whether $5 is [ACTG],X - i.e. two alternate alleles, the second 'X' - if yes it changes: (1) $5 to everything before ',' (gsub(/,.*/,"",$5) and (2) all genotype "2" to "." (i.e. missing)
	#than it checks whether $5 is X,[ACTG] - i.e. two alternate alleles, the first 'X' - if yes it changes: (1) $5 to everything after ',' (gsub(/.*,/,"",$5), (2) all genotypes "1" to "." (i.e. missing) and (3) all genotypes "2" to "1"
	#between every command there is sed command replacing ' ' to TABs
	#FIRST - this treats 'X' in $5 (i.e. ALT)
	awk -F"\t" '{ if ($5=="X") {$5="*";for (i=10; i<=NF; i++) if ($i=="1") {$i="."}; print $0} else {print $0} }' dataToModif.txt | sed "s/ /\t/g" | awk -F"\t" 'match($5, /[ACTG],X/) {gsub(/,.*/,"",$5);for (i=10; i<=NF; i++) if ($i=="2") {$i="."}; print $0; next}1' | sed "s/ /\t/g" | awk -F"\t" 'match($5, /X,[ACTG]/) {gsub(/.*,/,"",$5);for (i=10; i<=NF; i++) if ($i=="1") {$i="."}; for (i=10; i<=NF; i++) if ($i=="2") {$i="1"}; print $0; next}1' | sed "s/ /\t/g" > FIRSTdataModif.txt
	#SECOND - this treats 'X' in $4 (i.e. REF)
	#check if there is 'X' in $$ and two alternatives in $5 - then it put first alternative to 'a' and second alternative to 'b' and rename $4 and $5 by 'a' and 'b', respectively
	#then it checks the values in columns from 10 until the end (for (i=10; i<=NF; i++)) - i.e. genotype values - and change all "0"s to "." (i.e. missing), all "1"s to "0"s and all "2"s to "1"s
	#next awk command checks if there is 'X' in $4 and names $4 by '*'
	awk -F"\t" '{if ($5~/[ACTG],[ACTG]/ && $4=="X") {a=(sprintf(substr($5,1,1)));b=(sprintf(substr($5,3,1))); $4=a; $5=b; for (i=10; i<=NF; i++) if ($i=="0") {$i="."}; for (i=10; i<=NF; i++) if ($i=="1") {$i="0"}; for (i=10; i<=NF; i++) if ($i=="2") {$i="1"}; print $0; } else {print $0} }' FIRSTdataModif.txt | sed "s/ /\t/g" | awk -F"\t" '{ if ($4=="X") {$4="*";for (i=10; i<=NF; i++) if ($i=="0") {$i="."}; print $0} else {print $0} }' | sed "s/ /\t/g" > SECONDdataModif.txt
	#remove unvariable sites (with '*', i.e. sites where the only variability is missing site)
	grep -v "*" SECONDdataModif.txt > FINALdataModif.txt
	#combine with header
	cat header.txt FINALdataModif.txt > concatenatedModifFinal.vcf
	rm header.txt data.txt dataToModif.txt FIRSTdataModif.txt SECONDdataModif.txt FINALdataModif.txt
	#filter the final VCF using BCFtools
	gzip concatenatedModifFinal.vcf
	#divide 'snpmiss' by '100' and filter with bcftools
	fltr=`echo "scale=2; $snpmiss / 100" | bc | sed 's/^\./0./'`
	bcftools view -e "AC==0 || AC==AN || F_MISSING > $fltr" -m2 -M2 -O z -o FILTEREDconcatenated_${snpmiss}missing.vcf.gz concatenatedModifFinal.vcf.gz
fi

#CONTINUE HERE !!!
#Create 0/1 matrix
echo -e "Creating 0/1 matrix...\n"
#unzip VCF
gzip -d -c FILTEREDconcatenated_${snpmiss}missing.vcf.gz > FILTEREDconcatenated_${snpmiss}missing.vcf
#only take columns 10+ (i.e., genotype data only)
grep -v "^##" FILTEREDconcatenated_${snpmiss}missing.vcf | cut -f2,10- > SNPmatrix.txt
#transpose the matrix (lines will be samples)
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' SNPmatrix.txt > SNPmatrixTransposed.txt

#Plot heatmap in R
if [[ $PBS_O_HOST == *".cz" ]]; then
	module unload python/3.7.7-gcc-8.3.0-t4loj4a #to avoid conflict with R4.4 (requires python3.9.12)
	module add r/4.4.0-gcc-10.2.1-ssuwpvb
	export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"
fi
echo -en "Plotting heatmap in R..."
if [[ $snpmiss -eq 0 ]]; then
	echo -e "using DistBinary(PERMANOVA) - no missing data\n"
	R --slave -f SNPheatmap.R 0 >> R.log 2>&1
else
	echo -e "using custom R function - missing data allowed\n"
	R --slave -f SNPheatmap.R 1 >> R.log 2>&1
fi

if [ ! -z "$selection" ]; then
	mv SNPheatmap.pdf ${selection}_SNPheatmap_${snpmiss}missing.pdf
	mv distmat.txt ${selection}_SNPdistmat_${snpmiss}missing.txt
else
	mv SNPheatmap.pdf ${data}_SNPheatmap_${snpmiss}missing.pdf
	mv distmat.txt ${data}_SNPdistmat_${snpmiss}missing.txt
fi

#Copy results to home
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNPheatmap/${snpmiss}missing
	cp *.{pdf,gz,fasta,txt,log} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNPheatmap/${snpmiss}missing
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNPheatmap/${snpmiss}missing
	cp *.{pdf,gz,fasta,txt,log} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNPheatmap/${snpmiss}missing
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/SNPheatmap/${snpmiss}missing
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/SNPheatmap/${snpmiss}missing
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08n
fi

echo -e "\nHybPhyloMaker8n finished...\n"

