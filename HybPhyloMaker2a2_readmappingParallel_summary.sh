#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker2a2_readmappingParallel_summary
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker2_readmappingParallel_summary
#$ -o HybPhyloMaker2_readmappingParallel_summary.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                 Script 02a2 - Read mapping in parallel summary               *
# *                                   v.1.8.0c                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2025 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker2a2 is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add datamash-1.3 #for data summary
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker2a2 is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir02a2
	cd workdir02a2
	#Add necessary modules
	module load datamash #not yet on Hydra
else
	echo -e "\nHybPhyloMaker2a2 is running locally...\n"
	echo -e "This summary is only for cluster environment. Exiting...\n"
	exit 3
fi

#Write log
logname=HPM2a2
echo -e "HybPhyloMaker2a2: summary of read mapping in parallel" > ${logname}.log
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
for set in data cp mapping probes cpDNACDS cpDNA mappingmethod conscall mincov majthres plurality nrns; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -en "Working with cpDNA"
	type="cp"
	cp $source/$cpDNACDS .
	probes=$cpDNACDS
elif [[ $cp =~ "full" ]]; then
	echo -e "Working with full plastomes\n"
	type="fullplastomeT"
elif [[ $cp =~ "no" ]]; then
	echo -en "Working with exons"
	type="exons"
	cp $source/$probes .
else
	echo -e "Variable 'cp' is not set to one of allowed values (yes, not, full). Exiting...\n"
	exit 3
fi
echo -n

#Copy list of samples
cp $path/10rawreads/SamplesFileNames.txt .
#Copy all fasta files
cp $path/$type/21mapped_${mappingmethod}/*.fasta .
#Copy all mapping summaries
cp $path/$type/21mapped_${mappingmethod}/mapping_summary*.txt .
if [[ $type =~ "exon" || $type =~ "cp" ]]; then
	#Copy all per target coverages
	mkdir coverage
	cp $path/$type/21mapped_${mappingmethod}/coverage/*_perTarget.txt coverage/
fi
#Rename the probe file and move it (just for the case it has the suffix *.fasta which would make problems in the next step)
if [[ $cp =~ "yes" ]]; then
	mv $cpDNACDS coverage/probes.fa
elif [[ $cp =~ "no" ]]; then
	mv $probes coverage/probes.fa
fi

#Creating a summary table
echo -e "\nCreating summary tables..."
#Produce tab-separated table
#Write headers (number of reads, nr paired reads, nr forward unpaired reads, nr reverse unpaired reads, nr mapped reads, % of mapped reads)
echo -e "Sample no.\tGenus\tSpecies\tTotal nr. reads\tNr. paired reads\tNr. forward unpaired reads\tNr. reverse unpaired reads\tNr. mapped reads\tPercentage of mapped reads" > tmpmap
cat mapping_summary*.txt >> tmpmap
mv tmpmap mapping_summary_${type}.txt
cp mapping_summary_${type}.txt $path/$type/21mapped_${mappingmethod}/

#Combine all fasta files into one
echo -e "\nCombining consensus FASTA files..."
cat *.fasta > consensus.fasta
#Remove line breaks from fasta file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' consensus.fasta > tmp && mv tmp consensus.fasta
mkdir -p $path/$type/30consensus
if [[ $cp =~ "yes" ]]; then
	cp consensus.fasta consensus_cpDNA.fasta
	cp consensus_cpDNA.fasta $path/$type/30consensus
elif [[ $cp =~ "full" ]]; then
	cp consensus.fasta consensus_fullplastome.fasta
	cp consensus_fullplastome.fasta $path/$type/30consensus
else
	cp consensus.fasta $path/$type/30consensus
fi

#Calculate percentage of missing data per accession (for fullplastome consensus)
if [[ $cp =~ "full" ]]; then
	#Calculate length of alignment: 1. get second line and count length, 2. decrease value by one (because previous command also counted LF)
	length=$(cat consensus.fasta | head -n 2 | tail -n 1 | wc -c)
	length=`expr $length - 1`
	#Replace newline with ' ' if line starts with '>' (i.e., merge headers with data into single line separated by space)
	cat consensus.fasta | sed '/^>/{N; s/\n/ /;}' > consensus.modif.fasta
	#Cut first part until space, i.e. header, and remove '>'
	cat consensus.modif.fasta | cut -f1 -d" " | sed 's/>//' > headers.txt
	#Cut only part after the first space, i.e., only sequence, change all missing data (-, ?, N) to 'n', replace all other characters then 'n' by nothing and print percentage of 'n's in each sequence
	cat consensus.modif.fasta | cut -f2 -d" " | sed 's/[?N-]/n/g' | sed 's/[^n]//g' | awk -v val=$length '{ print (length*100)/val }' > missingpercentage.txt
	paste headers.txt missingpercentage.txt > fullplastomes_missingperc.txt
	#Calculate mean of all values
	echo -e "MEAN\t$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' fullplastomes_missingperc.txt)" > mean.txt
	cat fullplastomes_missingperc.txt mean.txt > tmp && mv tmp fullplastomes_missingperc.txt
	rm consensus.modif.fasta headers.txt missingpercentage.txt mean.txt
	#Copy results to home
	cp fullplastomes_missingperc.txt $path/$type/30consensus
fi

#Calculate ambiguous bases (number and percentage) if ConsensusFixer was used
if [[ $conscall =~ "consensusfixer" ]]; then
	echo -e "\nCalculating ambiguous bases for ConsensusFixer results..."
	#print length of sequences (i.e., print line length if the first character on line is not '>')
	awk '{ if (substr($1,1,1) !~ /^>/ ) print length($0) }' consensus.fasta > totallength.txt
	#Remove '?'s (introduced when N in reference), only in sequences (i.e., not on lines starting with '>')
	sed -i '/>/!s/\?//g' consensus.fasta
	#print length of sequences without '?'
	awk '{ if (substr($1,1,1) !~ /^>/ ) print length($0) }' consensus.fasta > totallengthNoQ.txt
	#Remove 'N's (introduced when too low coverage), only in sequences (i.e., not on lines starting with '>')
	sed -i '/>/!s/N//g' consensus.fasta
	#Replace newline with ' ' if line starts with '>' (i.e., merge headers with data into single line separated by space)
	cat consensus.fasta | sed '/^>/{N; s/\n/ /;}' > consensus.modif.fasta
	#Cut first part until space, i.e. header, and remove '>'
	cat consensus.modif.fasta | cut -f1 -d" " | sed 's/>//' > headers.txt
	#Calculate sequence length for each line
	cat consensus.modif.fasta | cut -f2 -d" " | awk '{ print length}' > lengths.txt
	#Cut only part after the first space, i.e., only sequence, modify to big letters only, change all ambiguities (RYSWKMBDHV) to 'x', replace all other characters then 'x' by nothing
	cat consensus.modif.fasta | cut -f2 -d" " | tr [a-z] [A-Z] | sed 's/[RYSWKMBDHV]/x/g' | sed 's/[^x]//g' > x.txt
	#Combine lengths.txt (sequence lengths) and x.txt (number of ambiguities) and print number and percentage of 'x's (i.e., ambiguities) in each sequence
	paste totallength.txt totallengthNoQ.txt lengths.txt x.txt | awk '{ print $1 "\t" $2 "\t" $3 "\t" length($4) "\t" (length($4)*100)/$3 }' > ambigpercentage.txt
	
	#Calculate number and percentage of each ambiguous base
	for base in R Y S W K M B D H V; do
		echo -e "Nr. $base\tPerc $base" > ${base}_percentage.txt
		cat consensus.modif.fasta | cut -f2 -d" " | sed "s/$base/x/g" | sed 's/[^x]//g' > ${base}x.txt
		paste lengths.txt ${base}x.txt | awk '{ print length($2) "\t" (length($2)*100)/$1 }' >> ${base}_percentage.txt
	done
	sed '1 i sample' headers.txt > headers2.txt
	paste headers2.txt *_percentage.txt > ambigbaseperc.txt
	paste headers.txt ambigpercentage.txt > ambigperc.txt
	#Calculate mean of all values
	echo -e "Sample\tTotalLength\tLengthWithout?\tCalledBases\tNrAmbig\tPercAmbig" > first.txt
	echo -e "MEAN\t$(awk '{ sum += $2; sum2 += $3; sum3 += $4; sum4 += $5; sum5 += $6; n++ } END { if (n > 0) print sum / n "\t" sum2 / n "\t" sum3 / n "\t" sum4 / n "\t" sum5 / n}' ambigperc.txt)" > mean.txt
	cat first.txt ambigperc.txt mean.txt > tmp && mv tmp ambigperc.txt
	rm consensus.modif.fasta headers.txt headers2.txt ambigpercentage.txt mean.txt *_percentage.txt lengths.txt *x.txt totallength.txt totallengthNoQ.txt first.txt
	#extract text after last '/' in $data (whole $data in no '/')
	data1=$(echo $data | rev | cut -d"/" -f1 | rev)
	cp ambigperc.txt $path/$type/30consensus/${data1}_ambigperc.txt
	cp ambigbaseperc.txt $path/$type/30consensus/${data1}_ambigbaseperc.txt
fi

# Make summary table for coverage
if [[ $type =~ "exon" || $type =~ "cp" ]]; then
	cd coverage
	echo -e "\nCreating summary tables for coverage..."
	ls *perTarget.txt | cut -d"." -f1 | sed 's/_perTarget//' > ListOfCoverageFiles.txt
	cat `ls *perTarget.txt | head -n 1` | cut -f 2,3,4 > infocolumn.txt
	echo -e "taxon\n5\n10\n20\n30\n50\n100" >> firstcolumn.txt
	for i in $(cat ListOfCoverageFiles.txt)
	do
		# Print header (species name)
		echo "$i" >> ${i}_exon_meancoverage.txt
		# Extract 7th column containing mean coverage per exon
		cat ${i}_perTarget.txt | cut -f7 | sed '1d' >> ${i}_exon_meancoverage.txt
		# Write samples name to file with 
		echo $i > ${i}_nrexons.txt
		for j in 5 10 20 30 50 100
		do
			cat ${i}_perTarget.txt | cut -f7 | awk -v val=$j ' NR>1 {if ($1>val) print $1;}' | wc -l >> ${i}_nrexons.txt
		done
	done
	# Combine information from all samples together
	paste infocolumn.txt *_exon_meancoverage.txt > exon_meancoverageALL.txt
	paste firstcolumn.txt *_nrexons.txt > numberexonsALL.txt
	# Copy summary table to home
	cp exon_meancoverageALL.txt $path/$type/21mapped_${mappingmethod}/coverage
	cp numberexonsALL.txt $path/$type/21mapped_${mappingmethod}/coverage
	# Calculate mean values per locus (a mean from per-exon means)
	# Delete first three rows (i.e., leaving only mean values)
	cut -f4- exon_meancoverageALL.txt | datamash transpose > coverage.txt
	# Prepare header
	echo locus | tr "\n" "\t" > header.txt #print 'locus' and replaces EOL by TAB, i.e., next line prints on the same line
	grep ">" probes.fa | cut -d'_' -f2 | datamash transpose >> header.txt
	echo -e "exon\t" | perl -0pe 's/\n\Z//' >> header.txt #print 'exonTAB' and removes very last character which is EOL, i.e., next line prints on the same line
	grep ">" probes.fa | cut -d'_' -f4 | datamash transpose >> header.txt
	cat header.txt coverage.txt > exon_meancoverageALLmodif.txt
	rm header.txt coverage.txt
	# Calculate mean per locus from per-exon values
	lines=$(wc -l < "exon_meancoverageALLmodif.txt")
	datamash -W transpose < "exon_meancoverageALLmodif.txt" | datamash -H groupby 1 mean 3-"$lines" | datamash transpose > locus_meancoverageALL.txt
	#modify output
	sed -i.bak 's/GroupBy(//' locus_meancoverageALL.txt
	sed -i.bak2 's/mean(//' locus_meancoverageALL.txt
	sed -i.bak3 's/)//' locus_meancoverageALL.txt
	rm *bak*
	cp exon_meancoverageALLmodif.txt $path/$type/21mapped_${mappingmethod}/coverage
	cp locus_meancoverageALL.txt $path/$type/21mapped_${mappingmethod}/coverage
	cd ..
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/21mapped_${mappingmethod}/

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir02a2
fi

echo -e "\nScript HybPhyloMaker2a2 finished...\n"

