#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=12:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=80gb
#PBS -j oe
#PBS -N HybPhyloMaker0d_create_test_dataset
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q mThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker0d_create_test_dataset
#$ -o HybPhyloMaker0d_create_test_dataset.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                       Script 0d - Preapare test dataset                      *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2023 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker0d is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	#. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add samtools-1.9
	module add bedtools
	module add seqtk
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker0d is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir0d
	cd workdir0d
	#Add necessary modules
	module load bioinformatics/samtools/1.3
else
	echo -e "\nHybPhyloMaker0d is running locally...\n"
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir0d
	cd workdir0d
fi

#Test if 'workdir' exist
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir0d)" ]; then
		echo -e "Directory 'workdir0d' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir0d/ 2>/dev/null
		exit 3
	fi
fi


#SET variables
#nrexons=250
echo -e "Nr. exons\t$nrexons"
echo
#probes=curcuma_HybSeqProbes.fa
#nrns=400

#copy data and probe reference
cp -r ${path}/exons/21mapped_${mappingmethod}/*.bam .
cp ${source}/${probes} .
sed -i 's/\x0D$//' $probes #remove 'CR' characters

#calculate other variables
nrex=$(($nrexons*2)) #number of lines containing nrexons is twice nrexons
#cut filename before first '.', i.e. remove suffix - does not work if there are other dots in reference file name
name=`ls $probes | cut -d'.' -f 1`
probes2=${name}_first${nrexons}.fa
name2=${name}_first${nrexons}

#create reference of desired number of exons
echo "Creating pseudoreferences..."
echo
head -n ${nrex} ${probes} > ${name}_first${nrexons}.fa

#create full pseudoreference (to get its length)
#print header to fasta file
echo ">${name}_with${nrns}Ns_beginend" > ${name}_with${nrns}Ns_beginend.fas
#print N $nrns times to variable $a
a=$(printf "%0.sN" $(seq 1 $nrns))
# 1. awk command to remove EOLs from lines not beginning with '>', i.e. all lines containing sequence for particular record are merged, i.e. each record on two lines only
# 2. sed command to delete all lines starting with '>', i.e. only sequences remain
# 3. tr command to replace all EOLs ("\n") by Q
# 4. sed command to replace all Qs by a sequence of $nrns Ns (saved in variable $a)
# 5. awk command to print $nrns Ns to the beginning and the end of the reference
cat $probes | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed '/>/d' | tr "\n" "Q" | sed "s/Q/$a/g" | awk -v val=$a '{ print val $0 val }' >> ${name}_with${nrns}Ns_beginend.fas
#full pseudoreference length (get line not containing '>' and count number of characters
tl=$(grep -v ">" ${name}_with${nrns}Ns_beginend.fas | wc -c)

#create pseudoreference from nrexons only
#print header to fasta file
echo ">${name2}_with${nrns}Ns_beginend" > ${name2}_with${nrns}Ns_beginend.fas
#a=$(printf "%0.sN" $(seq 1 $nrns))
cat $probes2 | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed '/>/d' | tr "\n" "Q" | sed "s/Q/$a/g" | awk -v val=$a '{ print val $0 val }' >> ${name2}_with${nrns}Ns_beginend.fas
#pseudoreference length for first nrexons
ll=$(grep -v ">" ${name2}_with${nrns}Ns_beginend.fas | wc -c)
#pseudoreference length for unused part (i.e. length after the last used exons: tl-ll)
ul=`expr $tl - $ll + 1`
#copy reduced probes and pseudoreference to home
cp $probes2 $source
echo -e "Probe file named '$probes2' was prepared and saved to '$source'..."
echo
cp ${name2}_with${nrns}Ns_beginend.fas $source
echo -e "Pseudoreference named '${name2}_with${nrns}Ns_beginend.fas' was prepared and saved to '$source'..."
echo


#create BED file for first nrexons
echo "Creating BED file..."
header=${name}_with${nrns}Ns_beginend
#delete lines not starting with '>'
sed '/^>/ d' $probes2 > reflines.txt
#delete empty lines
sed -i.bak '/^$/d' reflines.txt
#count length of each line
awk '{print length($0)}' reflines.txt > lengths.txt
#set $linesum to $nrns + 1 (i.e., starting position of the first exon)
linesum=`expr $nrns + 1`
cat lengths.txt | while read line
do
	#print BED file (name and two numbers: 1. start ($linesum) 2. end of the exon ($linesum + length of the exon - 1)
	echo $header $linesum `expr $linesum + $line - 1` >> ${name2}_with${nrns}Ns_beginend_exon_positions.bed
	#increaese $linesum by length of the region and $nrns (so now it is start position of the next exon)
	linesum=`expr $linesum + $line + $nrns`
done
#replace spaces by tabs
cat ${name2}_with${nrns}Ns_beginend_exon_positions.bed | tr ' ' '\t' > tmp && mv tmp ${name2}_with${nrns}Ns_beginend_exon_positions.bed
#remove some files
rm reflines* lengths.txt
echo

#copy BED file to home
cp ${name2}_with${nrns}Ns_beginend_exon_positions.bed $source
echo -e "\nBED file '${name2}_with${nrns}Ns_beginend_exon_positions.bed' was prepared and saved to '$source'..."

#make dir for results
mkdir 10rawreads_first${nrexons}

#loop over BAM samples
ls *.bam > samples.txt #create lost of BAM files
for file in $(cat samples.txt); do
	filename=`ls $file | cut -d'.' -f 1`
	echo Reducing ${filename}
	#extract reads for first nrexons from full BAM file using BED file of nrexons
	samtools view -b -L ${name2}_with${nrns}Ns_beginend_exon_positions.bed $file > ${filename}_first${nrexons}.bam
	#number of reads mapped to first nrexons
	nrused=$(samtools view -c ${filename}_first${nrexons}.bam)
	#extract reads after the last used exon, i.e. $ul-$tl (=unused reads)
	samtools sort -o ${filename}_sorted.bam $file
	samtools index ${filename}_sorted.bam
	samtools view -b ${filename}_sorted.bam "${header}:${ll}-${tl}" > ${filename}_unused.bam
	#number of reads apped to regions after the last used exon
	nrunused=$(samtools view -c ${filename}_unused.bam)
	#sort BAM files by name (necessary for later extraction of paired FASTQ files)
	samtools sort -n -o ${filename}_first${nrexons}.nsort.bam ${filename}_first${nrexons}.bam
	samtools sort -n -o ${filename}_unused.nsort.bam ${filename}_unused.bam
	#export reads as FASTQ (only full pairs are exported, single reads discarded)
	bedtools bamtofastq -i ${filename}_first${nrexons}.nsort.bam -fq ${filename}_first${nrexons}.R1.fastq -fq2 ${filename}_first${nrexons}.R2.fastq 2>/dev/null
	bedtools bamtofastq -i ${filename}_unused.nsort.bam -fq ${filename}_unused.R1.fastq -fq2 ${filename}_unused.R2.fastq 2>/dev/null
	#exported nr reads in first nrexons (i.e., nr reads in FASTQ file)
	nrused2=$((`wc -l < ${filename}_first${nrexons}.R1.fastq` / 4))
	echo -e "\t" ${nrused2} "mapped reads"
	#subsample unused reads (get a half of used)
	unusedtoadd=`expr $nrused2 / 2`
	seqtk sample -s 123 ${filename}_unused.R1.fastq $unusedtoadd > ${filename}_unusedSub.R1.fastq
	seqtk sample -s 123 ${filename}_unused.R2.fastq $unusedtoadd > ${filename}_unusedSub.R2.fastq
	echo -e "\t" $unusedtoadd "unmapped reads"
	#combine selected and subsamples unused reads
	cat ${filename}_first${nrexons}.R1.fastq ${filename}_unusedSub.R1.fastq > ${filename}_L001_R1_001_first${nrexons}.fastq
	cat ${filename}_first${nrexons}.R2.fastq ${filename}_unusedSub.R2.fastq > ${filename}_L001_R2_001_first${nrexons}.fastq
	#gzip fastq files
	gzip ${filename}_L001_R?_001_first${nrexons}.fastq
	#make dir for results
	mkdir 10rawreads_first${nrexons}/${filename}
	#move resulted FASTQ files
	mv ${filename}_L001_R?_001_first${nrexons}.fastq.gz 10rawreads_first${nrexons}/${filename}
	#delete unused files (everything starting with $filename, results already moved)
	if [[ ! ${filename} == "" ]]; then
		rm ${filename}*
	fi
done

#copy all results back
cp -r 10rawreads_first${nrexons} ${path}

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir0d
fi

echo -e "\nScript HybPhyloMaker0d finished...\n"
