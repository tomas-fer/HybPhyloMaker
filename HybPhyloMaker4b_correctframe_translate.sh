#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=4:mem=4gb:scratch_local=4gb
#PBS -j oe
#PBS -N HybPhyloMaker4b_translate
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 4
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker4b_translate
#$ -o HybPhyloMaker4b_translate.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *       Script 04b - Put exons to correct reading frame and translate          *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************


#Take MAFFT exon alignment files, check for correct reading frame (+ remove incomplete triplets at the beginning and end of each exon)
#and translate to amino acids using EMBOSS
#Produce concatenation of
#(1) corrected exonic nucleotide alignments (+ generate 2 partition files - by exons and by codon positions within exons)
#(2) translated exons alignments

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker4b is running on MetaCentrum..."
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add emboss-6.5.7
	module add perl-5.10.1
	#module add python-3.4.1-gcc
	#module add debian8-compat
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker4b is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir -p workdir04b
	cd workdir04b
	#Add necessary modules
	module load bioinformatics/emboss/6.6.0
	module load bioinformatics/anaconda3/2.3.0 #for python3
else
	echo -e "\nHybPhyloMaker4b is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir -p workdir04b
	cd workdir04b
fi

#Write log
logname=HPM4b
echo -e "HybPhyloMaker4b: put exons to correct reading frame and translate" > ${logname}.log
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
for set in data cp otherpslx othersource maxstop; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
else
	echo -e "Working with exons\n"
	type="exons"
fi

#Copy scripts
cp $source/catfasta2phyml.pl .
cp $source/AMAS.py .
#Copy MAFFT exon alignments
cp $path/$type/60mafft/*.mafft .
#Make a dir for corrected exons and concatenated loci
mkdir $path/$type/61mafft_corrected
mkdir $path/$type/62mafft_translated
mkdir $path/$type/80concatenated_exon_alignments_corrected
mkdir $path/$type/90concatenated_exon_alignments_translated

#Translate into frame 1, 2 and 3, count number of stop codons and save it to file
echo -ne "Translating exons to all three reading frames and selecting the correct one..."
ls *.mafft > listOfMAFFTFiles.txt
for mafftfile in $(cat listOfMAFFTFiles.txt)
do
	#Replace gaps in sequences by 'n'
	sed -i.bak '/>/!s/-/n/g' $mafftfile
	echo $mafftfile >> ${mafftfile}_stopcodonnr_overview.txt
	for i in 1 2 3
	do
		transeq fasta::$mafftfile -frame $i stdout 2>/dev/null | fgrep -o \* | wc -l >> stopcodonnr_${mafftfile}.txt
	done
	cat stopcodonnr_${mafftfile}.txt >> ${mafftfile}_stopcodonnr_overview.txt
	cp ${mafftfile}_stopcodonnr_overview.txt $path/$type/62mafft_translated
	#Select the frame with least number of stop codons
	min=$(sort -n stopcodonnr_${mafftfile}.txt | head -1)
	frame=$(grep -n ^${min}$ stopcodonnr_${mafftfile}.txt | cut -c1 | head -1)
	#Generate modified DNA alignment (introduce frame shift)
	#Number of added character is $frame-1 (i.e., frame1=0, frame2=1, frame3=2)
	if [ "$frame" = 2 ]; then
		add=2
	elif [ "$frame" = 3 ]; then
		add=1
	elif [ "$frame" = 1 ]; then
		add=0
	fi
	addcharacter=$(printf '%*s' $add '' | tr ' ' 'n')
	#Add 'n's at the beginning of each nucleotide sequence (first remove line breaks from FASTA file - awk command, second add '???' after each line not beginning with '>' - sed command)
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${mafftfile} |  sed '/>/!s/^/'"$addcharacter"'/' > ${mafftfile}.corrframe
	#Add '???' at the end of each nucleotide sequence (i.e., after each line not beginning with '>')
	#sed -i '/>/!s/$/???/' ${mafftfile}.corrframe
	#Remove first three characters from sequences in case of frame shift to remove incomplete triplet
	if [ "$frame" != 1 ]; then
		sed -i.bak '/>/!s/^...//' ${mafftfile}.corrframe
	fi
	#Remove incomplete triplets from the end of exon (0, 1 or 2 bases)
	#Compute length of the alignment (length of the first accession, i.e., second line)
	len=$(cat ${mafftfile}.corrframe | head -2 | tail -1 | awk '{ print length }')
	remove=`expr $len % 3` #Calculate modulo after division by 3, i.e., number of bases to be remove`
	#Remove last 0, 1 or 2 characters from sequences (i.e., from each line not beginning with '>')
	sed -i.bak '/>/!s/.\{'"$remove"'\}$//' ${mafftfile}.corrframe
	#Copy corrected exon alignment to home
	cp ${mafftfile}.corrframe $path/$type/61mafft_corrected
	
	#Translate corrected exon
	transeq fasta::$mafftfile.corrframe -frame 1 -outseq ${mafftfile}.translated 2>/dev/null
	#remove line breaks from translated FASTA file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${mafftfile}.translated > tmp && mv tmp ${mafftfile}.translated
	#Add '???' at the end of each protein sequence (add '???' after each line not beginning with '>' - sed command)
	#sed -i.bak '/>/!s/$/???/' ${mafftfile}.translated
	#Remove last two characters from fasta headers (i.e., each line beginning with '>') - were introduced by EMBOSS transeq
	sed -i.bak '/>/s/.\{2\}$//' ${mafftfile}.translated
	#Copy translated exon alignment to home
	cp ${mafftfile}.translated $path/$type/62mafft_translated
done
echo -e "done\n"

echo -en "Preparing summary table..."
#combine the information about stop codons by frame
#paste *overview.txt > stop_codons_by_frame.txt #this usually does not work with more loci than 4096
touch stop_codons_by_frame.txt
for f in *overview.txt; do
        cat stop_codons_by_frame.txt | paste - $f > temp
        cp temp stop_codons_by_frame.txt
done
cut -f2- stop_codons_by_frame.txt > tmp && mv tmp stop_codons_by_frame.txt
rm temp


#transpose the matrix
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
}' stop_codons_by_frame.txt > stop_codons_by_frame_transposed.txt
awk '{print $1 }' stop_codons_by_frame_transposed.txt > names.txt
awk '{print $2 " " $3 " " $4 }' stop_codons_by_frame_transposed.txt > stop_codons_by_frame.txt
#put the lowest number as last column
awk '{min=$1;for(i=1;i<=NF;i++){if($i < min) min = $i} print $0 " " min}' stop_codons_by_frame.txt > tmp && mv tmp stop_codons_by_frame.txt
#put the number of zero stop codons as last column
awk '{zeros=0;for(i=1;i<=NF-1;i++){if($i=="0") zeros++ } print $0 " " zeros }' stop_codons_by_frame.txt > tmp && mv tmp stop_codons_by_frame.txt
#combine back with names
paste names.txt stop_codons_by_frame.txt > tmp && mv tmp stop_codons_by_frame.txt
cat stop_codons_by_frame.txt | tr " " "\t" > tmp && mv tmp stop_codons_by_frame.txt
echo -e "done\n"

#Filter exons
echo -ne "Filtering exons..."
#1. Remove exons with more than one reading frame with 0 stop codons
awk '{if($6>1) print $1}' stop_codons_by_frame.txt > removed_more_than_1_possible_reading_frame.txt
awk '{if($6<=1) print $0}' stop_codons_by_frame.txt > selected_exons.txt
#2. Remove exons with lowest number of stop codons exceeding threshold (maxstop)
awk -v val=$maxstop '{if($5>val) print $1}' stop_codons_by_frame.txt > removed_lowest_number_of_stop_codons_exceeded_maxstop.txt
awk -v val=$maxstop '{if($5<=val) print $1}' selected_exons.txt > tmp && mv tmp selected_exons.txt
#Copy lists to home
cp selected_exons.txt $path/$type/62mafft_translated
cp removed*.txt $path/$type/62mafft_translated

#Modify translation summary and copy to home
sed -i.bak 's/To_align_Assembly_//g' stop_codons_by_frame.txt
sed -i.bak2 's/.fasta.mafft//g' stop_codons_by_frame.txt
cp stop_codons_by_frame.txt $path/$type/62mafft_translated

#Copy selected exons to folder
mkdir selected
mkdir selectedtranslated
for i in $(cat selected_exons.txt); do
	cp ${i}.translated selectedtranslated
	cp ${i}.corrframe selected
done

#Modify selected exons (replace stopcodons by NNN)
for i in $(cat selected_exons.txt); do
	sed 's/.../&+/g;s/+$//' selected/${i}.corrframe | sed 's/taa/NNN/g' | sed 's/tag/NNN/g' | sed 's/tga/NNN/g' | sed 's/+//g' > selected/tmp && mv selected/tmp selected/${i}.corrframe
	cp selected/${i}.corrframe $path/$type/61mafft_corrected
done
echo -e "done\n"
#-----------------------CONCATENATE MODIFIED NUCLEOTIDE EXON ALIGNMENTS-----------------------
echo -ne "Concatenating corrected exons..."
#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft to To_align_Assembly_10372_*mafft)
#(all files starting with "To_align_Assembly_10372_" will be merged)
ls -1 selected/*.corrframe | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft.corrframe/g > fileNamesForConcat.txt
#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. Assembly_10372
ls -1 selected/*.corrframe | cut -d'_' -f4 | sort -u | sed s/^/CorrectedAssembly_/g > fileNamesForSaving.txt
#Combine both files (make single file with two columns)
paste fileNamesForConcat.txt fileNamesForSaving.txt > fileForLoop.txt

#Concatenate the exon alignments (values from first and second column of fileForLoop.txt are assigned to variable 'a' and 'b', respectively),
#transform fasta to phylip format, copy results from scratch to home
cat fileForLoop.txt | while read -r a b; do
	python3 AMAS.py concat -i selected/$a -f fasta -d dna -u fasta -t selected/${b}.fasta >/dev/null
	python3 AMAS.py concat -i selected/$a -f fasta -d dna -u phylip -t selected/${b}.phylip >/dev/null
	#Modify partition file for partitioning by codon positions
	# cat partitions.txt			#take the file
	# cut -d"_" -f2- 				#take 2nd and following parts of each line separated by "_" (remove 'p1_' etc. introduced by AMAS)
	# awk '{ print $0 "\\3;" }'		#
	# awk '{for(i=0;i<3;i++)print}'	#print each line 3-times
	# awk -F '[=-]' '{ if (NR % 3 == 1) print $1 "_1=" $2 "-" $3; else if (NR % 3 == 2) print $1 "_2=" $2+1 "-" $3; else print $1 "_3=" $2+2 "-" $3}'
	#								#set separators to '=' and '-' and on 1st, 2nd and 3rd line do different modifications (add _1, _2 or _3; leave unchanged, add 1 or add 2 to $2, i.e. start position)
	# sed 's/-/ - /g'				#replace "-" by " - "
	# sed 's/=/ = /g'				#replace "=" by " = "
	cat partitions.txt | cut -d"_" -f2- | awk '{ print $0 "\\3;" }' | awk '{for(i=0;i<3;i++)print}' | awk -F '[=-]' '{ if (NR % 3 == 1) print $1 "_1=" $2 "-" $3; else if (NR % 3 == 2) print $1 "_2=" $2+1 "-" $3; else print $1 "_3=" $2+2 "-" $3}' | sed 's/-/ - /g' | sed 's/=/ = /g' | sed -e 's/^/DNA, /g' -e 's/To_align_//g' | sed 's/;//g' > selected/${b}.codonpart.file
	#modify and rename partitions.txt
	cat partitions.txt | cut -d"_" -f2- | sed -e 's/^/DNA, /g' -e 's/To_align_//g' > selected/${b}.part
	#perl catfasta2phyml.pl -f $a > $b.fasta
	#perl catfasta2phyml.pl $b.fasta > $b.phylip
	#Copy concatenated alignments and partition files to home
	cp selected/$b.* $path/$type/80concatenated_exon_alignments_corrected
done
echo -e "done\n"

#-----------------------CONCATENATE TRANSLATED EXON ALIGNMENTS-----------------------
echo -ne "Concatenating translated exons..."
#Modify mafft file names (from, i.e., To_align_Assembly_10372_Contig_1_516.fasta.mafft.translated to To_align_Assembly_10372_*mafft.translated)
#(all files starting with "To_align_Assembly_10372_" will be merged)
ls -1 selectedtranslated/*.translated | cut -d'_' -f4 | sort -u | sed s/^/To_align_Assembly_/g | sed s/\$/_*mafft.translated/g > fileNamesForConcat.txt
#Modify mafft file names - prepare names for saving concatenate alignments (not possible to use names above as they contain '*'), e.g. Assembly_10372
ls -1 selectedtranslated/*.translated | cut -d'_' -f4 | sort -u | sed s/^/TranslatedAssembly_/g > fileNamesForSaving.txt
#Combine both files (make single file with two columns)
paste fileNamesForConcat.txt fileNamesForSaving.txt > fileForLoop.txt
#Concatenate the exon alignments (values from first and second column of fileForLoop.txt are assigned to variable 'a' and 'b', respectively),
#transform fasta to phylip format, copy results from scratch to home
cat fileForLoop.txt | while read -r a b; do
	python3 AMAS.py concat -i selectedtranslated/$a -f fasta -d aa -u fasta -t selectedtranslated/${b}.fasta >/dev/null
	python3 AMAS.py concat -i selectedtranslated/$a -f fasta -d aa -u phylip -t selectedtranslated/${b}.phylip >/dev/null
	#modify and rename partitions.txt
	sed -i.bak -e 's/^/WAG, /g' -e 's/To_align_//g' partitions.txt
	mv partitions.txt selectedtranslated/${b}.part
	#perl catfasta2phyml.pl -f $a > $b.fasta
	#perl catfasta2phyml.pl $b.fasta > $b.phylip
	#Copy translated concatenated alignments to home
	cp selectedtranslated/$b.* $path/$type/90concatenated_exon_alignments_translated
done
echo -e "done\n"

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/80concatenated_exon_alignments_corrected

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir04b
fi

echo -e "\nScript HybPhyloMaker04b finished...\n"