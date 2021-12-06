#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N HybPhyloMaker2a_readmappingParallel
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker2_readmappingParallel
#$ -o HybPhyloMaker2_readmappingParallel.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                     Script 02a - Read mapping in parallel                    *
# *                                   v.1.8.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker2a is running on MetaCentrum...\n"
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add perl-5.10.1
	module add gcc-4.8.4
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker2a is running on Hydra...\n"
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir02a
	cd workdir02a
	#Add necessary modules
else
	echo -e "\nHybPhyloMaker2a is running locally...\n"
	echo -e "...parallel processing is not supported, run HybPhyloMaker1_rawprocess.sh instead"
	echo -e "Exiting..."
	exit 3
fi

#Test if 'workdir' exist
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir02a)" ]; then
		echo -e "Directory 'workdir02a' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir02a/ 2>/dev/null
		exit 3
	fi
fi

#Advice to run HybPhyloMaker2 if mapping=no and exit
if [[ $mapping =~ "no" ]]; then
	echo -e "Variable 'mapping'=no. Run HybPhyloMaker2_readmapping.sh for consensus call only. Exiting...\n"
	exit 3
fi

#Write log
logname=HPM2a
echo -e "HybPhyloMaker2a: mapping reads in parallel" > ${logname}.log
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
for set in data cp mapping probes cpDNACDS mappingmethod conscall mincov majthres plurality nrns; do
	printf "%-25s %s\n" `echo -e "${set}:\t" ${!set}` >> ${logname}.log
done

#Setting for the case when working with cpDNA
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
	cp $source/$cpDNACDS .
	
else
	echo -e "Working with exons\n"
	type="exons"
	cp $source/$probes .
fi

#Test if folders for results exist
if [ -d "$path/$type/21mapped_${mappingmethod}" ] && [[ $mapping =~ "yes" ]]; then
	echo -e "Directory '$path/$type/21mapped_${mappingmethod}' already exists. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir02a/ 2>/dev/null
	exit 3
elif [ -d "$path/$type/30consensus" ]; then
	echo -e "Directory '$path/$type/30consensus' already exists. Delete it or rename before running this script again. Exiting...\n"
	rm -d ../workdir02a/ 2>/dev/null
	exit 3
fi

#Test data structure
echo -en "Testing input data structure..."
if [ -f "$path/10rawreads/SamplesFileNames.txt" ]; then
	#Copy SamplesFileNames.txt and modify it
	cp $path/10rawreads/SamplesFileNames.txt .
	#Add LF at the end of last line in SamplesFileNames.txt if missing
	sed -i.bak '$a\' SamplesFileNames.txt
	#Delete empty lines from SamplesFileNames.txt (if any)
	sed -i.bak2 '/^$/d' SamplesFileNames.txt
	rm *bak*
	for sample in $(cat SamplesFileNames.txt); do
		if [ $mapping == "yes" ]; then
			if [ ! -d "$path/20filtered/$sample" ]; then #Test if each samples-specific folder exists
				echo -e "Directory $sample does not exist within '20filtered'.\n"
				rm SamplesFileNames.txt
				rm -d ../workdir02a/ 2>/dev/null
				exit 3
			else
				if [ ! -f "$path/20filtered/$sample/${sample}-1P_no-dups.fastq" ] || [ ! -f "$path/20filtered/$sample/${sample}-2P_no-dups.fastq" ] || [ ! -f "$path/20filtered/$sample/${sample}-1U" ] || [ ! -f "$path/20filtered/$sample/${sample}-2U" ]; then #Test if filtered FASTQ files exist
					if [ ! -f "$path/20filtered/$sample/${sample}-1P_no-dups.fastq.gz" ] || [ ! -f "$path/20filtered/$sample/${sample}-2P_no-dups.fastq.gz" ] || [ ! -f "$path/20filtered/$sample/${sample}-1U.gz" ] || [ ! -f "$path/20filtered/$sample/${sample}-2U.gz" ]; then #Test if filtered and compressed FASTQ files exist
						echo -e "Appropriate filtered fastq files missing in $sample folder...\n"
						rm SamplesFileNames.txt
						rm -d ../workdir02a/ 2>/dev/null
						exit 3
					else
						compressed=yes
					fi
				fi
			fi
		else
			if [ ! -f "$path/$type/21mapped_${mappingmethod}/${sample}.bam" ]; then
				echo -e "$sample.bam does not exist within '$type/21mapped_${mappingmethod}'.\n"
				rm SamplesFileNames.txt
				rm -d ../workdir02a/ 2>/dev/null
				exit 3
			fi
		fi
	done
	if [[ $cp =~ "yes" ]]; then
		if [ ! -f "$source/$cpDNACDS" ]; then
			echo -e "'$cpDNACDS' is missing in 'HybSeqSource'. Exiting...\n"
			rm -d ../workdir02a/ 2>/dev/null
			exit 3
		else
			name=$(echo $cpDNACDS | cut -d"." -f1)
			if [ ! -f "$source/${name}_with${nrns}Ns_beginend.fas" ]; then
				echo -e "${name}_with${nrns}Ns_beginend.fas does not exist within 'HybSeqSource'. Run HybPhyloMaker0b_preparereference.sh first.\n"
				rm SamplesFileNames.txt
				rm -d ../workdir02a/ 2>/dev/null
				exit 3
			fi
		fi
	else
		if [ ! -f "$source/$probes" ]; then
			echo -e "$probes does not exist within 'HybSeqSource'.\n"
			rm SamplesFileNames.txt
			rm -d ../workdir02a/ 2>/dev/null
			exit 3
		else
			name=$(echo $probes | cut -d"." -f1)
			if [ ! -f "$source/${name}_with${nrns}Ns_beginend.fas" ]; then
				echo -e "${name}_with${nrns}Ns_beginend.fas does not exist within 'HybSeqSource'. Run HybPhyloMaker0b_preparereference.sh first.\n"
				rm SamplesFileNames.txt
				rm -d ../workdir02a/ 2>/dev/null
				exit 3
			fi
		fi
	fi
else
	echo -e "List of samples (SamplesFileNames.txt) is missing. Should be in 10rawreads...\n"
	rm -d ../workdir02a/ 2>/dev/null
	exit 3
fi
echo -e "OK for running HybPhyloMaker2\n\n"

#Make a new folder for results
mkdir -p $path/$type
if [ ! -d "$path/$type/21mapped_${mappingmethod}" ]; then
	mkdir $path/$type/21mapped_${mappingmethod}
fi
if [ ! -d "$path/$type/21mapped_${mappingmethod}/coverage" ]; then
	mkdir $path/$type/21mapped_${mappingmethod}/coverage
fi

#Create new 'submitRawProcessJobs.sh' and make it executable
touch $path/$type/21mapped_${mappingmethod}/submitMappingJobs.sh
chmod +x $path/$type/21mapped_${mappingmethod}/submitMappingJobs.sh

#Copy list of samples
cp $path/10rawreads/SamplesFileNames.txt .

for file in $(cat SamplesFileNames.txt); do
	echo -e "Processing sample $file..."
	echo '#!/bin/bash' >> ${file}.sh
	echo '#----------------MetaCentrum----------------' >> ${file}.sh
	echo '#PBS -l walltime=2:0:0' >> ${file}.sh
	echo '#PBS -l select=1:ncpus=4:mem=4gb:scratch_local=16gb' >> ${file}.sh
	echo '#PBS -j oe' >> ${file}.sh
	echo '#PBS -N readmapping_for_'"${file}" >> ${file}.sh
	echo '#PBS -m abe' >> ${file}.sh
	echo '#-------------------HYDRA-------------------' >> ${file}.sh
	echo '#$ -S /bin/bash' >> ${file}.sh
	echo '#$ -pe mthread 4' >> ${file}.sh
	echo '#$ -q mThC.q' >> ${file}.sh
	echo '#$ -l mres=4G,h_data=4G,h_vmem=4G' >> ${file}.sh
	echo '#$ -cwd' >> ${file}.sh
	echo '#$ -j y' >> ${file}.sh
	echo '#$ -N readmapping_for_'"${file}" >> ${file}.sh
	echo '#$ -o readmapping_for_'"${file}"'.log' >> ${file}.sh
	echo '#Complete path and set configuration for selected location' >> ${file}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${file}.sh
	echo '  #Add necessary modules' >> ${file}.sh
	echo '  module add bowtie2-2.2.4' >> ${file}.sh
	echo '  module add bwa-0.7.15' >> ${file}.sh
	echo '  module add samtools-1.9' >> ${file}.sh
	echo '  module add perl-5.10.1' >> ${file}.sh
	echo '  module add gcc-4.8.4' >> ${file}.sh
	echo '  module add python34-modules-gcc #adds also kindel' >> ${file}.sh
	echo '  module add ococo-2016-11' >> ${file}.sh
	echo '  module add jdk-8' >> ${file}.sh
	echo '  cd $SCRATCHDIR' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  #Add necessary modules' >> ${file}.sh
	echo '  module load bioinformatics/bowtie2/2.2.9' >> ${file}.sh
	echo '  module load bioinformatics/bwa/0.7.12' >> ${file}.sh
	echo '  module load bioinformatics/samtools/1.3' >> ${file}.sh
	echo '  module load bioinformatics/anaconda3/5.1 #adds also kindel' >> ${file}.sh
	echo '  module load bioinformatics/fastuniq/1.1' >> ${file}.sh
	echo '  #module load bioinformatics/ococo/ #???' >> ${file}.sh
	echo '  module load java/1.8' >> ${file}.sh
	echo '  mkdir -p workdir02a_${file}' >> ${file}.sh
	echo '  cd workdir02a_${file}' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo 'path='"$path" >> ${file}.sh
	echo 'source='"$source" >> ${file}.sh
	echo 'type='"$type" >> ${file}.sh
	echo 'compressed='"$compressed" >> ${file}.sh
	echo 'file='"$file" >> ${file}.sh
	echo 'location='"$location" >> ${file}.sh
	echo 'data='"$data" >> ${file}.sh
	echo 'cp='"$cp" >> ${file}.sh
	echo 'mapping='"$mapping" >> ${file}.sh
	echo 'probes='"$probes" >> ${file}.sh
	echo 'cpDNACDS='"$cpDNACDS" >> ${file}.sh
	echo 'mappingmethod='"$mappingmethod" >> ${file}.sh
	echo 'conscall='"$conscall" >> ${file}.sh
	echo 'mincov='"$mincov" >> ${file}.sh
	echo 'majthres='"$majthres" >> ${file}.sh
	echo 'plurality='"$plurality" >> ${file}.sh
	echo 'nrns='"$nrns" >> ${file}.sh
	echo '#Copy programs' >> ${file}.sh
	echo 'if [[ $conscall =~ "consensusfixer" ]]; then' >> ${file}.sh
	echo '  cp $source/ConsensusFixer.jar .' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo '#Copy pseudoreference' >> ${file}.sh
	echo 'if [[ $cp =~ "yes" ]]; then' >> ${file}.sh
	echo '  probes=$cpDNACDS' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo 'name=$(echo $probes | cut -d"." -f1)' >> ${file}.sh
	echo 'cp $source/${name}_with${nrns}Ns_beginend.fas .' >> ${file}.sh
	echo '#Make index from pseudoreference' >> ${file}.sh
	echo 'if [[ $mappingmethod =~ "bowtie2" ]]; then' >> ${file}.sh
	echo '  echo -en "Indexing pseudoreference for bowtie2..."' >> ${file}.sh
	echo '  bowtie2-build ${name}_with${nrns}Ns_beginend.fas pseudoreference.index 1>indexing_pseudoreference.log' >> ${file}.sh
	echo '  cp indexing_pseudoreference.log $path/$type/21mapped_bowtie2/' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  echo -en "Indexing pseudoreference for bwa...\n"' >> ${file}.sh
	echo '  bwa index ${name}_with${nrns}Ns_beginend.fas 2>indexing_pseudoreference.log' >> ${file}.sh
	echo '  cp indexing_pseudoreference.log $path/$type/21mapped_bwa/' >> ${file}.sh
	echo '  fi' >> ${file}.sh
	echo '#Write headers for a summary table (number of reads, nr paired reads, nr forward unpaired reads, nr reverse unpaired reads, nr mapped reads, % of mapped reads)' >> ${file}.sh
	echo 'echo -e "Sample no.\tGenus\tSpecies\tTotal nr. reads\tNr. paired reads\tNr. forward unpaired reads\tNr. reverse unpaired reads\tNr. mapped reads\tPercentage of mapped reads" > mapping_summary_${file}.txt' >> ${file}.sh
	echo '#Read mapping' >> ${file}.sh
	echo 'if [[ $compressed =~ "yes" ]]; then' >> ${file}.sh
	echo '  suffix=".gz"' >> ${file}.sh
	echo '  #copy fastq.gz files and count number of reads' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-1P_no-dups.fastq.gz .' >> ${file}.sh
	echo '  nr1P=$(gunzip -c ${file}-1P_no-dups.fastq.gz | awk '\'{s++}END{print s/4}\'')' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-2P_no-dups.fastq.gz .' >> ${file}.sh
	echo '  nr2P=$(gunzip -c ${file}-2P_no-dups.fastq.gz | awk '\'{s++}END{print s/4}\'')' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-1U.gz .' >> ${file}.sh
	echo '  nr1U=$(gunzip -c ${file}-1U.gz | awk '\'{s++}END{print s/4}\'')' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-2U.gz .' >> ${file}.sh
	echo '  nr2U=$(gunzip -c ${file}-2U.gz | awk '\'{s++}END{print s/4}\'')' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  suffix=""' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-1P_no-dups.fastq .' >> ${file}.sh
	echo '  nr1P=$(awk '\'{s++}END{print s/4}\'' ${file}-1P_no-dups.fastq)' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-2P_no-dups.fastq .' >> ${file}.sh
	echo '  nr2P=$(awk '\'{s++}END{print s/4}\'' ${file}-2P_no-dups.fastq)' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-1U .' >> ${file}.sh
	echo '  nr1U=$(awk '\'{s++}END{print s/4}\'' ${file}-1U)' >> ${file}.sh
	echo '  cp $path/20filtered/${file}/${file}-2U .' >> ${file}.sh
	echo '  nr2U=$(awk '\'{s++}END{print s/4}\'' ${file}-2U)' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo 'if [[ $mappingmethod =~ "bowtie2" ]]; then' >> ${file}.sh
	echo '  echo "Mapping using bowtie2..."' >> ${file}.sh
	echo '  #Bowtie2 parameters are derived from --very-sensitive-local' >> ${file}.sh
	echo '  #set parameters for mapping using bowtie2' >> ${file}.sh
	echo '  score=G,20,8' >> ${file}.sh
	echo '  bowtie2 --local -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --score-min $score -x pseudoreference.index -1 ${file}-1P_no-dups.fastq${suffix}  -2 ${file}-2P_no-dups.fastq${suffix} -U ${file}-1U${suffix},${file}-2U${suffix} -S ${file}.sam 2>${file}_bowtie2_out.txt' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  echo "Mapping pair-end reads using bwa..."' >> ${file}.sh
	echo '  bwa mem ${name}_with${nrns}Ns_beginend.fas ${file}-1P_no-dups.fastq${suffix} ${file}-2P_no-dups.fastq${suffix} > ${file}_paired.sam 2>${file}_bwa_out.txt' >> ${file}.sh
	echo '  echo "Mapping orphaned reads using bwa..."' >> ${file}.sh
	echo '  if [[ $compressed =~ "yes" ]]; then' >> ${file}.sh
	echo '    gunzip ${file}-1U.gz' >> ${file}.sh
	echo '    gunzip ${file}-2U.gz' >> ${file}.sh
	echo '  fi' >> ${file}.sh
	echo '  cat ${file}-1U ${file}-2U > ${file}-unpaired' >> ${file}.sh
	echo '  bwa mem ${name}_with${nrns}Ns_beginend.fas ${file}-unpaired > ${file}_unpaired.sam 2>>${file}_bwa_out.txt' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo '#create BAM from SAM' >> ${file}.sh
	echo 'echo "Converting to BAM..."' >> ${file}.sh
	echo 'if [[ $mappingmethod =~ "bowtie2" ]]; then' >> ${file}.sh
	echo '  samtools view -bS -o ${file}.bam ${file}.sam 2>/dev/null' >> ${file}.sh
	echo '  rm ${file}.sam' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  samtools view -bS -o ${file}_paired.bam ${file}_paired.sam 2>/dev/null' >> ${file}.sh
	echo '  samtools view -bS -o ${file}_unpaired.bam ${file}_unpaired.sam 2>/dev/null' >> ${file}.sh
	echo '  samtools merge ${file}.bam ${file}_paired.bam ${file}_unpaired.bam 2>/dev/null' >> ${file}.sh
	echo '  rm ${file}_paired.sam ${file}_unpaired.sam ${file}_paired.bam ${file}_unpaired.bam' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo '#number of mapped reads in BAM' >> ${file}.sh
	echo 'nrmapped=$(samtools view -F 0x04 -c ${file}.bam)' >> ${file}.sh
	echo '#number of all reads in BAM' >> ${file}.sh
	echo 'nrall=$(samtools view -c ${file}.bam)' >> ${file}.sh
	echo '#calculate percentage of mapped reads' >> ${file}.sh
	echo 'percmapped=$(echo -e "scale=5;100 * ($nrmapped / $nrall)" | bc)' >> ${file}.sh
	echo '#extracting header data' >> ${file}.sh
	echo 'number=$(cut -d'\'_\'' -f2 <<< $file)' >> ${file}.sh
	echo 'genus=$(cut -d'\'-\'' -f1 <<< $file)' >> ${file}.sh
	echo 'species=$(cut -d'\'_\'' -f1 <<< $file | cut -d'\'-\'' -f2)' >> ${file}.sh
	echo '#adding data to table' >> ${file}.sh
	echo 'echo -e "$number\t$genus\t$species\t$nrall\t$nr1P\t$nr1U\t$nr2U\t$nrmapped\t$percmapped" > mapping_summary_${file}.txt' >> ${file}.sh
	echo '#sort and index' >> ${file}.sh
	echo 'echo "Sorting and indexing BAM..."' >> ${file}.sh
	echo 'samtools sort ${file}.bam -o ${file}_sorted.bam' >> ${file}.sh
	echo 'rm ${file}.bam' >> ${file}.sh
	echo 'mv ${file}_sorted.bam ${file}.bam' >> ${file}.sh
	echo 'samtools index ${file}.bam' >> ${file}.sh
	echo '#copy results to home' >> ${file}.sh
	echo 'cp ${file}.bam $path/$type/21mapped_${mappingmethod}' >> ${file}.sh
	echo 'cp ${file}.bam.bai $path/$type/21mapped_${mappingmethod}' >> ${file}.sh
	echo 'cp ${file}_${mappingmethod}_out.txt $path/$type/21mapped_${mappingmethod}' >> ${file}.sh
	echo 'cp mapping_summary_${file}.txt $path/$type/21mapped_${mappingmethod}' >> ${file}.sh
	echo '#CONSENSUS USING KINDEL/OCOCO/ConsensusFixer' >> ${file}.sh
	echo 'if [[ $conscall =~ "ococo" ]]; then' >> ${file}.sh
	echo '  echo "Making consensus with OCOCO..."' >> ${file}.sh
	echo '  ococo -i ${file}.bam -x ococo64 -c $mincov -F ${file}.fasta 2>/dev/null' >> ${file}.sh
	echo 'elif [[ $conscall =~ "consensusfixer" ]]; then' >> ${file}.sh
	echo '  echo "Making consensus with ConsensusFixer..."' >> ${file}.sh
	echo '  if [ ! -f "${file}.bam.bai" ]; then' >> ${file}.sh
	echo '    samtools sort ${file}.bam -o ${file}_sorted.bam' >> ${file}.sh
	echo '    rm ${file}.bam' >> ${file}.sh
	echo '    mv ${file}_sorted.bam ${file}.bam' >> ${file}.sh
	echo '    samtools index ${file}.bam' >> ${file}.sh
	echo '  fi' >> ${file}.sh
	echo '  #call consensus with ambiguous bases with ConsensusFixer' >> ${file}.sh
	echo '  java -jar ConsensusFixer.jar -i ${file}.bam -r ${name}_with${nrns}Ns_beginend.fas -plurality $plurality -mcc $mincov -dash' >> ${file}.sh
	echo '  #add EOL at the end of the file' >> ${file}.sh
	echo "  sed -i '\$a\' consensus.fasta" >> ${file}.sh
	echo "  #change '-' (introduce by ConsensusFixer when coverage is low) by 'N'" >> ${file}.sh
	echo "  sed -i 's/-/N/g' consensus.fasta" >> ${file}.sh
	echo '  mv consensus.fasta ${file}.fasta' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  echo "Making consensus with kindel..."' >> ${file}.sh
	echo '  kindel -m $mincov -t $majthres ${file}.bam > ${file}.fasta' >> ${file}.sh
	echo 'fi' >> ${file}.sh
	echo '#change name in fasta file' >> ${file}.sh
	echo "sed -i.bak '1d' \${file}.fasta #delete first line" >> ${file}.sh
	echo 'echo ">$file" > header.txt' >> ${file}.sh
	echo 'cat header.txt ${file}.fasta > tmp && mv tmp ${file}.fasta' >> ${file}.sh
	echo 'rm header.txt' >> ${file}.sh
	echo '#Remove line breaks from fasta file' >> ${file}.sh
	echo awk \''!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'\'' ${file}.fasta > tmp && mv tmp ${file}.fasta' >> ${file}.sh
	echo '#put $nrns Ns to variable '\'a\'' and $nrns ?s to variable '\'b\''' >> ${file}.sh
	echo 'a=$(printf "%0.sN" $(seq 1 $nrns))' >> ${file}.sh
	echo 'b=$(printf "%0.s?" $(seq 1 $nrns))' >> ${file}.sh
	echo '#replace all Ns separating exons by '?'' >> ${file}.sh
	echo 'sed -i.bak "s/$a/$b/g" ${file}.fasta' >> ${file}.sh
	echo '#copy results to home' >> ${file}.sh
	echo 'cp ${file}.fasta $path/$type/21mapped_${mappingmethod}' >> ${file}.sh
	echo '#Per exon coverage calculation' >> ${file}.sh
	echo 'echo -e "\nPer exon coverage calculation...\n"' >> ${file}.sh
	echo '#Copy references and BED file' >> ${file}.sh
	echo 'reference=${name}_with${nrns}Ns_beginend.fas' >> ${file}.sh
	echo 'cp $source/$reference .' >> ${file}.sh
	echo 'header=$(grep ">" ${name}_with${nrns}Ns_beginend.fas | sed "s/>//")' >> ${file}.sh
	echo 'bedfile=${header}_exon_positions.bed' >> ${file}.sh
	echo 'cp $source/$bedfile .' >> ${file}.sh
	echo '#Get picard tools' >> ${file}.sh
	echo 'echo -e "Downloading picard tools...\n"' >> ${file}.sh
	echo 'wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar 2> /dev/null' >> ${file}.sh
	echo '#Index reference' >> ${file}.sh
	echo 'samtools faidx ${reference}' >> ${file}.sh
	echo '#Create sequence dictionary from reference' >> ${file}.sh
	echo 'echo -e "Creating sequence dictionary from reference\n" > picard_${file}.log' >> ${file}.sh
	echo 'java -jar picard.jar CreateSequenceDictionary -R ${reference} -O reference.dict 2>> picard_${file}.log' >> ${file}.sh
	echo 'echo >> picard_${file}.log' >> ${file}.sh
	echo '#Convert BED file to IntervalList' >> ${file}.sh
	echo 'echo -e "Converting BED file to IntervalList\n" >> picard_${file}.log' >> ${file}.sh
	echo 'java -jar picard.jar BedToIntervalList -I ${bedfile} -O list -SD reference.dict 2>> picard_${file}.log' >> ${file}.sh
	echo 'echo >> picard_${file}.log' >> ${file}.sh
	echo 'echo -e "Calculating per exon coverage..."' >> ${file}.sh
	echo '#compute metrics (%GC, coverage)' >> ${file}.sh
	echo 'java -jar picard.jar CollectHsMetrics -I ${file}.bam -O ${file}_metrics.txt -R ${reference} -BAIT_INTERVALS list -TARGET_INTERVALS list -PER_TARGET_COVERAGE ${file}_perTarget.txt 2>> picard_${file}.log' >> ${file}.sh
	echo '# Copy results from SCRATCHDIR to HOME' >> ${file}.sh
	echo 'cp ${file}_perTarget.txt $path/$type/21mapped_${mappingmethod}/coverage' >> ${file}.sh
	echo 'cp picard_${file}.log $path/$type/21mapped_${mappingmethod}/coverage' >> ${file}.sh
	echo '#delete BAM' >> ${file}.sh
	echo 'rm ${file}.bam' >> ${file}.sh
	echo 'rm ${file}.bam.bai' >> ${file}.sh
	echo '#Clean scratch/work directory' >> ${file}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${file}.sh
	echo '  #delete scratch' >> ${file}.sh
	echo '  rm -rf $SCRATCHDIR/*' >> ${file}.sh
	echo 'else' >> ${file}.sh
	echo '  cd ..' >> ${file}.sh
	echo '  rm -r workdir02a_'"${file}" >> ${file}.sh
	echo 'fi' >> ${file}.sh
	
	chmod +x ${file}.sh
	if [[ $location == "1" ]]; then
		cp ${file}.sh $path/$type/21mapped_${mappingmethod}
		#qsub ${file}.sh
		echo 'qsub '"${file}"'.sh' >> $path/$type/21mapped_${mappingmethod}/submitMappingJobs.sh
	else
		cp ${file}.sh $path/$type/21mapped_${mappingmethod}
		cp ${file}.sh ..
		echo 'qsub '"${group}"'.sh' >> ../submitMappingJobs.sh
	fi
done

echo -e "done\n"

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
cp ${logname}.log $path/$type/21mapped_${mappingmethod}

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir02a
fi

echo -e "\nHybPhyloMaker 2a finished..."
if [[ $location == "2" ]]; then
	echo -e "\nGo to homedir and run submitMappingJobs.sh...\n"
else
	echo -e "\nGo to $path/$type/21mapped_${mappingmethod} and run submitMappingJobs.sh..."
	echo -e "This starts parallel mapping of reads to reference."
	echo -e "\nAfter all jobs finish run script HybPhyloMaker2a2 in order to summarize the mapped data...\n"
fi
