#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N cpRefFromGB
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThM.q
#$ -l mres=4G,h_data=4G,h_vmem=4G,himem
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8a2_Astral4
#$ -o HybPhyloMaker8a2_Astral4.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *            Script 0f - Prepare plastome reference from GenBank file          *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2024 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#optionally DOWNLOAD full plastome accession from GenBank (requires Entrez Direct)
#(https://www.ncbi.nlm.nih.gov/books/NBK179288/)

#EXTRACT sequences (in FASTA format) from full plastome provided in GenBank format
#- all features (CDS, tRNA, rRNA), separated to exons
#- all sequences among them (introns, spacers)

#PREPARE reference for HybPhyloMaker
#- modify names
#- delete regions shorter than $short (defined below!)

#KNOWN BUGS:
#1. limited working in case the '/gene' descriptor is missing in CDS description (like, e.g., in MK460223)
#   in this case this script requires 'product.txt' for proper modification to gene name ('/gene' descriptor is added)
#   however, some '/product' descriptions matching appropriate gene names have to be added/modified to 'product.txt'
#   known problem: 'RNA polymerase beta subunit' or 'RNA polymerase beta' describes rpoB, rpoC1 and rpoC2 - manual check necessary!!!
#   known problem: some CDS have '/product' description as 'no product string in file'
#   also, if tRNAs and rRNAs are not annotated in GB file, some of the non-coding regions extracted by this script include them
#2. Does not properly work for plastomes of parazitic plants, e.g. HG514460
#   in this case only remaining functional genes are extracted (those annotated as CDS)
#   non-coding regions then contain non-coding parts as well as non-functional genes (i.e. genes with '/pseudo')
#   tRNAs are extracted both functional and pseudo

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker0f is running on MetaCentrum..."
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
	module add entrezdirect
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker0f is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir0f
	cd workdir0f
	#Add necessary modules
	module load entrezdirect???
else
	echo -e "\nHybPhyloMaker0f is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir0f
	cd workdir0f
fi

if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir0f)" ]; then
		echo -e "Directory 'workdir0f' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir0f 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM0f
echo -e "HybPhyloMaker0f: Prepare plastome reference from GenBank file" > ${logname}.log
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
for set in cpGBfile cpGBnr mincplength; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

#Add necessary files
cp $source/product.txt .
#remove all comments from product.txt
grep -v "^#" product.txt > tmp && mv tmp product.txt

#Make dir for results
mkdir -p $path/cp/00reference

#download GenBank file using efetch
if [[ "$cpGBnr" == "" && "$cpGBfile" == "" ]]; then
	echo -e "Name of the GenBank file (cpGBfile) OR GenBank AccessionNumber (cpGBnr) not provided.\nExiting...\n" && exit 3
elif [[ -z $cpGBfile ]]; then
	echo -en "Trying to download Accession: ${cpGBnr} ... "
	#test if esearch returns Counts=1
	estest=$(esearch -db nuccore -query "$cpGBnr" | grep "Count" | awk -F'>|<' '{ print $3 }')
	if [[ ${estest} != 1 ]]; then
		echo -e "not found.\nExiting...\n" && exit 3
	else
		efetch -db nuccore -id ${cpGBnr} -format genbank > ${cpGBnr}.gb
		${cpGBfile}=${cpGBnr}.gb
		echo -e "completed\n"
	fi
elif [[ ! -f "$source/${cpGBfile}" ]]; then
	echo -e "File: ${cpGBfile} was not found in HybSeqSource.\nExiting...\n" && exit 3
else
	cp $source/${cpGBfile} .
fi

#set variables
name=`echo ${cpGBfile} | cut -d'.' -f1`
org=$(grep ORGANISM ${cpGBfile} | awk '{ print $2"-"$3 }')
org2=$(grep ORGANISM ${cpGBfile} | awk '{ print $2" "$3 }')
fastah=$(echo ">${org}_${name}")

echo -e "\nWorking with '${org2} (${name})'..."

#get just the sequence
echo -e "\nExtracting sequence..."
#take everything after 'ORIGIN', remove numbers, remove spaces, remove '//', remove CRs, remove LFs
sed -e '1,/ORIGIN/d' ${cpGBfile} | sed 's/[0-9]*//g' | sed 's/ //g' | sed 's/\/\///' | sed 's/\x0D$//' | tr -d '\n' > ${name}.fasta

#Modify GB file (in order to put join feature information on a single line)
#delete everything between ',' and any number [0-9], i.e. putting the join description of some CDS on a single line
perl -0777 -pe 's/(,)\s*([0-9])/$1$2/g' ${cpGBfile} > ${cpGBfile}.modif #'perl -0777' = do not any splitting

# Modify GB file in case there are no '/gene' descriptors for CDS
if [ "$(grep "/gene" ${cpGBfile} | wc -l)" -lt 20 ]; then
	echo "Converting '/products' to '/genes'..."
	# modify long names of '/product' (the case of 'ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit')
	perl -0777 -pe 's/(oxygenase)\s*(large)/$1 $2/g' ${cpGBfile}.modif > tmp && mv tmp ${cpGBfile}.modif
	# change '/product' to gene name according to product.txt
	cat product.txt | while read -r a b; do
		c=$(echo "$b" | sed 's/\//\\\//g') #change '/' to '\/' before sed on the next line
		sed -i "s/\"$c\"/\"$a\"/" ${cpGBfile}.modif
	done
	# delete EOL on lines with 'note'
	sed -i '/note/ { N; s/\n// }' ${cpGBfile}.modif
	# remove unnecessary lines
	sed -i '/\/note/ d' ${cpGBfile}.modif
	sed -i '/\/trans_splicing/ d' ${cpGBfile}.modif
	sed -i '/\/codon_start/ d' ${cpGBfile}.modif
	sed -i '/\/transl/ d' ${cpGBfile}.modif
	sed -i '/locus_tag/ d' ${cpGBfile}.modif
	sed -i '/\/exception/ d' ${cpGBfile}.modif
	# change '/product' to '/gene'
	sed -i 's/\/product/\/gene/' ${cpGBfile}.modif
fi
rm product.txt

echo "Extracting features..."
#select lines containing CDS, tRNA and one following line (the line with gene name)
#CDS is surrounded by spaces ( CDS ) to avoid capturing the string from AA translation part
grep -A1 -E " CDS |tRNA " ${cpGBfile}.modif > resultCDS.txt
#add '--'' as last line (to mimic grep output after later file merging)
#select lines containing rRNA and two following lines (the line with gene name)
echo "--" >> resultCDS.txt
#remove spaces in gene names if any
sed -i '/gene/s/ //g' resultCDS.txt
#replace two or more spaces by ',', remove lines starting with ',/locus' or ',/gene' or '/old_locus_tag' and select lines containing rRNA and one following line
#(this is to remove other descriptors and get product' descriptor as the next line to rRNA)
sed 's/ \{2,\}/,/g' ${cpGBfile}.modif | sed '/^,\/locus/ d' | sed '/^,\/gene/ d' | sed '/^,\/old_locus_tag/ d' | grep -A1 -E ",rRNA" --no-group-separator | sed "s/,rRNA/--\n,rRNA/" | sed 1d > resultR.txt
rm ${cpGBfile}.modif
#grep -A2 -E "rRNA " ${cpGBfile} > resultR.txt
sed -i 's/ ribosomal RNA//' resultR.txt
#replace spaces by ','
sed -i 's/ \{1,\}/,/g' resultCDS.txt
sed -i 's/ \{1,\}/,/g' resultR.txt
#delete first ',' at each line
sed -i 's/,//' resultCDS.txt
sed -i 's/,//' resultR.txt
#delete first '/' at each line
sed -i 's/\///' resultCDS.txt
sed -i 's/\///' resultR.txt
#select correct name by deleting other lines
sed -i '/^locus/ d' resultR.txt
#merge both files
cat resultCDS.txt resultR.txt > result.txt
rm resultCDS.txt resultR.txt
#replace EOLs by ','
cat result.txt | tr "\n" "," > tmp && mv tmp result.txt
#delete "(" and ")"
sed -i 's/(//g' result.txt
sed -i 's/)//g' result.txt
#delete 'complement'
sed -i 's/complement//g' result.txt
#add ',' after 'join'
sed -i 's/join/join,/g' result.txt
#replace ',--,' by EOL
sed -i "s/,--,/\n/g" result.txt
#replace '..' by ','
sed -i 's/\.\./,/g' result.txt
#delete 'gene='
sed -i 's/gene=//' result.txt
#delete 'product='
sed -i 's/product=//' result.txt
#remove double quotes
sed -i 's/\"//g' result.txt
#delete last character on last line (if it is an excessive ',')
sed -i '$s/,$//' result.txt
#add newline at the end of the file
sed -i '$a\' result.txt
#replace '-' by 'x'
sed -i 's/-/x/' result.txt

echo "Treating joined features..."
#extract lines containing join
grep join result.txt > join.txt
#calculates number of regions to join
awk -F "," '{ print ( NF - 3 ) / 2 }' join.txt > nrRegions.txt
#change ',' to TABs
cat join.txt | tr "," "\t" > tmp && mv tmp join.txt
#merge two files
paste nrRegions.txt join.txt > nrRegionsJoin.txt

#loop over genes with more parts and separates parts to one line each
cat nrRegionsJoin.txt | while read a b c d e f g h i j; do
	if [ "$a" -eq 2 ]; then
		echo -e "$b,$d,$e,${h}_1" >> joinResult.txt
		echo -e "$b,$f,$g,${h}_2" >> joinResult.txt
	elif [ "$a" -eq 3 ]; then
		echo -e "$b,$d,$e,${j}_1" >> joinResult.txt
		echo -e "$b,$f,$g,${j}_2" >> joinResult.txt
		echo -e "$b,$h,$i,${j}_3" >> joinResult.txt
	fi
done
rm join.txt nrRegions.txt nrRegionsJoin.txt

#remove lines with 'join' from results.txt
sed -i '/join/d' result.txt
#combine result.txt and joinResult.txt
cat result.txt joinResult.txt > finalResult.txt
rm result.txt joinResult.txt
#change ',' to TABs
cat finalResult.txt | tr "," "\t" > tmp && mv tmp finalResult.txt
#sort according 2nd field (i.e., starting position)
sort -n -k2 finalResult.txt > finalResultSorted.txt
rm finalResult.txt

#add intergenic regions
echo "Adding intergenic regions..."
#generate file starting positions of intergenic regions
awk '{ print $3 + 1 "\t" $4}' finalResultSorted.txt > startIntergenic.txt
#get name of the last gene
last=$(tail -1 finalResultSorted.txt | cut -f4)
#add as first line
sed -i "1i 1\t$last" startIntergenic.txt
#merge sorted file and 
paste finalResultSorted.txt startIntergenic.txt > finalPlusIntergenic.txt
rm startIntergenic.txt
#write noncoding regions
awk '{ print "non\t" $5 "\t" $2 - 1 "\t" $6 "-" $4}' finalPlusIntergenic.txt > intergenic.txt
rm finalPlusIntergenic.txt
#get the length of the plastome
length=$(head -n1 ${cpGBfile} | sed 's/ \{1,\}/,/g' | cut -d',' -f3)
#get the end of last coding gene and increase by one (start position of last non-coding region)
lastcodpos=$(tail -n1 finalResultSorted.txt | cut -f3)
lastcodpos=$((lastcodpos+1))
#delete last line
head -n -1 intergenic.txt > tmp && mv tmp intergenic.txt
#create last line
echo -e "non\t$lastcodpos\t$length\t$last-end" > lastline.txt
cat intergenic.txt lastline.txt > tmp && mv tmp intergenic.txt
rm lastline.txt

#merge with coding and non-coding and sort
cat finalResultSorted.txt intergenic.txt | sort -n -k2 > tmp && mv tmp finalResultSorted.txt
rm intergenic.txt

#treat negative-length non-coding regions (established by the script when genes are overlapping)
echo "Treating negative-lengths regions..."
#calculate region lengths and add it as forth column
awk '{ print $1 "\t" $2 "\t" $3 "\t" $3 - $2 + 1 "\t" $4 }' finalResultSorted.txt > tmp && mv tmp finalResultSorted.txt
#print only lines with non-negative region length (4th column > 0)
awk '{ if ($4 > 0) print $0 }' finalResultSorted.txt > tmp && mv tmp finalResultSorted.txt

#create reference
echo "Creating reference..."
x=1
cat finalResultSorted.txt | while read a b c d e; do
	echo -e ">${x}_${x}_${e}_${a}" >> ${name}_reference.fasta
	cut -c $b-$c ${name}.fasta >> ${name}_reference.fasta
	x=$((x+1))
done

#Add organism name to fasta record
sed -i "1i ${fastah}" file${name}.fasta

echo "Removing empty records..."
#RS declares records separated by '>', if the record contains a second line $2 (i.e. is not an empty fasta record), print the record (add a > in front).
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ${name}_reference.fasta > tmp && mv tmp ${name}_reference.fasta

echo -e "feature\tstart\tstop\tlength\tname" > header.txt
cat header.txt finalResultSorted.txt > ${name}_plastomeFeaturePositions.txt
rm header.txt finalResultSorted.txt

#modifying reference for the use with HybPhyloMaker
echo "Modifying reference for HybPhyloMaker..."
#1. separate to types (tRNA, rRNA, CDS, non, CDS+non)
grep "tRNA$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_tRNA.fasta
grep "rRNA$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_rRNA.fasta
grep "CDS$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_CDS.fasta
grep "non$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_non.fasta
grep -e "CDS$" -e "non$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_CDSnon.fasta
#2. modify names
#first sed command replaces last occurrence of '_' by 'zzz'
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_tRNA.fasta > ${name}_reference_tRNAHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_rRNA.fasta > ${name}_reference_rRNAHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_CDS.fasta > ${name}_reference_CDSHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_non.fasta > ${name}_reference_nonHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_CDSnon.fasta > ${name}_reference_CDSnonHPM.fasta
#3. remove regions shorter than 200 bp
grep "^[^>].\{$mincplength\}" -B1 --no-group-separator ${name}_reference_tRNAHPM.fasta > ${name}_reference_tRNAHPM${mincplength}.fasta
grep "^[^>].\{$mincplength\}" -B1 --no-group-separator ${name}_reference_rRNAHPM.fasta > ${name}_reference_rRNAHPM${mincplength}.fasta
grep "^[^>].\{$mincplength\}" -B1 --no-group-separator ${name}_reference_CDSHPM.fasta > ${name}_reference_CDSHPM${mincplength}.fasta
grep "^[^>].\{$mincplength\}" -B1 --no-group-separator ${name}_reference_nonHPM.fasta > ${name}_reference_nonHPM${mincplength}.fasta
grep "^[^>].\{$mincplength\}" -B1 --no-group-separator ${name}_reference_CDSnonHPM.fasta > ${name}_reference_CDSnonHPM${mincplength}.fasta

#statistics
echo "Calculating nr of regions..."
#full reference
echo -e "CDS\nnon\ntRNA\nrRNA" > ${name}_statHeader.txt
grep "CDS$" ${name}_reference.fasta | wc -l > ${name}_statD.txt
grep "non$" ${name}_reference.fasta | wc -l >> ${name}_statD.txt
grep "tRNA$" ${name}_reference.fasta | wc -l >> ${name}_statD.txt
grep "rRNA$" ${name}_reference.fasta | wc -l >> ${name}_statD.txt
paste ${name}_statHeader.txt ${name}_statD.txt > ${name}_statT.txt
rm ${name}_statHeader.txt ${name}_statD.txt

#after short removal
grep "CDS$" ${name}_reference_CDSHPM${mincplength}.fasta | wc -l > ${name}_statD.txt
grep "non$" ${name}_reference_nonHPM${mincplength}.fasta | wc -l >> ${name}_statD.txt
grep "tRNA$" ${name}_reference_tRNAHPM${mincplength}.fasta | wc -l >> ${name}_statD.txt
grep "rRNA$" ${name}_reference_rRNAHPM${mincplength}.fasta | wc -l >> ${name}_statD.txt
paste ${name}_statT.txt ${name}_statD.txt > ${name}_statF.txt
rm ${name}_statT.txt ${name}_statD.txt
echo -e "type\tfull\tlonger${mincplength}" > header.txt
cat header.txt ${name}_statF.txt > ${name}_stat.txt
rm header.txt ${name}_statF.txt

#Copy results to home
cp *.{gb,fasta,txt} $path/cp/00reference
#Copy HPM references to HybSeqSource
cp *HPM*.fasta $source
cp ${cpGBfile} $source

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/cp/00reference

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir0f
fi

echo -e "\nScript HybPhyloMaker0f finished...\n"
