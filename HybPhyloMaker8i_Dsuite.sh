#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker8i_Dsuite
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -pe mthread 8
#$ -q mThC.q
#$ -l mres=3G,h_data=3G,h_vmem=3G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker8i_Dsuite
#$ -o HybPhyloMaker8i_Dsuite.log

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                             Script 08i - Dsuite                              *
# *                                   v.1.8.0d                                   *
# *                  Roman Ufimov, Martha Kandziora & Tomas Fer                  *
# *       Dept. of Botany, Charles University, Prague, Czech Republic, 2022      *
# *                           tomas.fer@natur.cuni.cz                            *
# ********************************************************************************


#Fast calculation of Paterson's D (ABBA-BABA) and the f4-ratio statistics across species
#according to https://github.com/millanek/Dsuite
#Malinsky M., Matschiner M & Svardal H (2020): Dsuite - Fast D-statistics and related admixture evidence from VCF files. Molecular Ecology Resources 21: 584-595.

#Take all gene alignments, concatenate them using AMAS and creates VCF file of variable sites using SNP-sites
#Uses ASTRAL species tree created by the script HybPhyloMaker8a_astral.sh

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker8i is running on MetaCentrum..."
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
	module add newick-utils-13042016
	module add conda-modules-py37 #to add SNP-sites later
	module add bcftools-1.11
	module add vcftools-0.1.16
	module unload python-2.6.6-gcc #avoid conflict with the next module
	#module add python36-modules-gcc
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker8i is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08i
	cd workdir08i
	#Add necessary modules
	module load bioinformatics/anaconda3/5.1 #adds NewickUtilities
	module load ??? #add SNP-sites on Hydra (maybe via conda?)
	module add bcftools???
	module add vcftools???
else
	echo -e "\nHybPhyloMaker8i is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir -p workdir08i
	cd workdir08i
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
echo -e "\nTesting if input data are available...not yet implemented"
#test for folder with alignments and for ASTRAL species tree
#Test if folder for results exits
if [[ $update =~ "yes" ]]; then
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Dsuite" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Dsuite' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08i 2>/dev/null
		exit 3
	fi
else
	if [ -d "${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Dsuite" ]; then
		echo -e "Directory '$path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Dsuite' already exists. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08i 2>/dev/null
		exit 3
	fi
fi
if [[ ! $location == "1" ]]; then
	if [ "$(ls -A ../workdir08i)" ]; then
		echo -e "Directory 'workdir08i' already exists and is not empty. Delete it or rename before running this script again. Exiting...\n"
		rm -d ../workdir08i 2>/dev/null
		exit 3
	fi
fi

#Write log
logname=HPM8i
echo -e "HybPhyloMaker8i: Dsuite" > ${logname}.log
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
for set in data selection cp corrected update tree OUTGROUP SNPs; do
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
wget https://raw.githubusercontent.com/mmatschiner/tutorials/master/analysis_of_introgression_with_snp_data/src/plot_d.rb 2>/dev/null
wget https://raw.githubusercontent.com/mmatschiner/tutorials/master/analysis_of_introgression_with_snp_data/src/plot_f4ratio.rb 2>/dev/null

# Make new dir for results
if [[ $update =~ "yes" ]]; then
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Dsuite
else
	mkdir $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Dsuite
fi

# Copy ASTRAL species tree
echo -e "Copying ASTRAL species tree...\n"
if [[ $update =~ "yes" ]]; then
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Astral/Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
else
	cp ${path}/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Astral/Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre .
fi

# Modify species tree
mv Astral_${MISSINGPERCENT}_${SPECIESPRESENCE}.tre sptree.tre
# replace ' ' back to '-' and '_'
sed -i 's/ \([^ ]*\) / \1_/g' sptree.tre #replace every second occurrence of ' ' by '_'
sed -i 's/ /-/g' sptree.tre #replace all spaces by '-'
# remove bootstrap values
nw_topology -bI sptree.tre > tmp && mv tmp sptree.tre
# rename ${OUTGROUP} to 'Outgroup' (requirement of Dsuite)
sed -i "s/$OUTGROUP/Outgroup/g" sptree.tre

# Create species set (each sample is a separate species)
nw_labels -I sptree.tre > SpeciesSet.txt
paste SpeciesSet.txt SpeciesSet.txt > tmp && mv tmp SpeciesSet.txt

# Concatenate gene alignments and create VCF
echo -e "Copying and modifying FASTA files...\n"
# copy all FASTA files (not selected only!)
find $path/${alnpath}/ -maxdepth 1 -name "*ssembly_*.fasta" -exec cp -t . {} + #to avoid 'Argument list too long' error
# concatenate FASTA files
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
# remove line breaks from FASTA
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' concatenated.fasta > tmp && mv tmp concatenated.fasta
sed -i 's/_contigs.fas//' concatenated.fasta
# rename ${OUTGROUP} to 'Outgroup' (requirement of Dsuite)
sed -i "s/$OUTGROUP/Outgroup/" concatenated.fasta
# put the outgroup as a first accession in FASTA
grep -A 1 "Outgroup" concatenated.fasta > final.fasta
# remove EOL after 'Outgroup', i.e. put 'Outgroup' and its sequence on a single line
perl -pe 's/Outgroup\n/Outgroup/g' concatenated.fasta > tmp && mv tmp concatenated.fasta
# take everything except the outgroup to final file
grep -v "Outgroup" concatenated.fasta >> final.fasta
mv final.fasta concatenated.fasta
# replace '?'s by 'N'
sed -i 's/\?/N/g' concatenated.fasta
# replace 'N's by '-'
sed -i '/^>/!s/N/-/g' concatenated.fasta

echo -e "Running SNP-site...\n"
if [[ $PBS_O_HOST == *".cz" ]]; then
	conda activate SNP-sites
	snp-sites -v concatenated.fasta > concatenated.vcf
	conda deactivate
	module unload conda-modules-py37
else
	snp-sites -v concatenated.fasta > concatenated.vcf
fi

#VCF parsing
grep ^# concatenated.vcf > header.txt
grep -v "^#" concatenated.vcf > snps.txt
cut -d'=' -f2 partitions.txt | tr '-' '\t' > bed.txt
#Calculate nr SNPs per exon
cat bed.txt | while read -r a b; do
	echo -en "from ${a} to ${b}\t" >> nrSNPsPerExonBefore.txt
	awk -v a=$a -v b=$b 'a<=$2 && $2<=b' snps.txt | wc -l >> nrSNPsPerExonBefore.txt
done

# Filter VCF
echo -e "Filtering VCF..."
gzip concatenated.vcf
bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2' -m2 -M2 -O z -o concatenated_bcf.vcf.gz concatenated.vcf.gz
gzip -d -c concatenated_bcf.vcf.gz > concatenated_bcf.vcf
grep ^# concatenated_bcf.vcf > headerF.txt
grep -v "^#" concatenated_bcf.vcf > snpsF.txt
rm concatenated_bcf.vcf
mv concatenated.vcf.gz FULLconcatenated.vcf.gz
mv concatenated_bcf.vcf.gz FILTEREDconcatenated.vcf.gz
#Calculate nr SNPs per exon after filtering
cat bed.txt | while read -r a b; do
	awk -v a=$a -v b=$b 'a<=$2 && $2<=b' snpsF.txt | wc -l >> nrSNPsPerExonFiltered.txt
done
paste nrSNPsPerExonBefore.txt nrSNPsPerExonFiltered.txt > nrSNPsPerExon.txt
echo -e "Interval\tBeforeFiltering\tAfterFiltering" > first.txt
cat first.txt nrSNPsPerExon.txt > tmp && mv tmp nrSNPsPerExon.txt
rm nrSNPsPerExonBefore.txt nrSNPsPerExonFiltered.txt first.txt

#VCF thinning
#first SNP per exon
cat bed.txt | while read -r a b; do
	awk -v a=$a -v b=$b 'a<=$2 && $2<=b' snpsF.txt | head -n 1 >> firstSNP.txt
done
cat headerF.txt firstSNP.txt > firstSNP.vcf
gzip firstSNP.vcf
#random SNP per exon
cat bed.txt | while read -r a b; do
	awk -v a=$a -v b=$b 'a<=$2 && $2<=b' snpsF.txt | sort -R > snpsF1.txt
	head -n 1 snpsF1.txt >> randomSNP.txt
done
rm snpsF1.txt
cat headerF.txt randomSNP.txt > randomSNP.vcf
gzip randomSNP.vcf
#thinning (i.e., one SNP in 100bp window)
vcftools --gzvcf FILTEREDconcatenated.vcf.gz --thin 100 --recode --out thinned
#nr thinned SNPs
grep -v "^#" thinned.recode.vcf > snps_thin.txt
mv thinned.recode.vcf thinned.vcf
gzip thinned.vcf
cat bed.txt | while read -r a b; do
	awk -v a=$a -v b=$b 'a<=$2 && $2<=b' snps_thin.txt | wc -l >> nrSNPsPerExonThin.txt
done
sed -i '1i Thinned' nrSNPsPerExonThin.txt
paste nrSNPsPerExon.txt nrSNPsPerExonThin.txt > tmp && mv tmp nrSNPsPerExon.txt

#SNP statistics
initial=$(wc -l < snps.txt)
filtered=$(wc -l < snpsF.txt)
onePerExon=$(wc -l < firstSNP.txt)
thinned=$(wc -l < snps_thin.txt)
echo -e "Initial\t${initial}\nFiltered\t${filtered}\nOnePerExon\t${onePerExon}\nThinned\t${thinned}" > SNPstat.txt
rm bed.txt snps.txt snpsF.txt header.txt headerF.txt firstSNP.txt randomSNP.txt snps_thin.txt nrSNPsPerExonThin.txt

#Prepare/select final VCF for Dsuite
if [[ $SNPs =~ "thinning" ]]; then
	cp thinned.vcf.gz Dsuite.vcf.gz
elif [[ $SNPs =~ "first" ]]; then
	cp firstSNP.vcf.gz Dsuite.vcf.gz
elif [[ $SNPs =~ "random" ]]; then
	cp randomSNP.vcf.gz Dsuite.vcf.gz
else
	echo "No valid option for SNP thinning provided. Exiting..." && exit 3
fi

# Run Dsuite
# compile
echo -e "\nRunning Dsuite..."
if [[ $PBS_O_HOST == *".cz" ]]; then
	#git clone https://github.com/millanek/Dsuite.git
	#IMPORTANT: if the analysis is not working with the most recent version of Dsuite, comment previous line and uncomment next three lines to get older version
	wget https://github.com/millanek/Dsuite/archive/60356e8493fa0bf1b85acea3be1b72df9dfb5881.zip
	unzip 60356e8493fa0bf1b85acea3be1b72df9dfb5881.zip
	mv Dsuite-* Dsuite
	cd Dsuite/
	make -j2
	cd ..
fi

#Set working python environment
module unload python/3.7.7-gcc-8.3.0-t4loj4a
module add python36-modules-gcc
# BBAA
echo -e "\nCreating BBAA..."
if [[ $PBS_O_HOST == *".cz" ]]; then
	Dsuite/Build/Dsuite Dtrios -c -n gene_flow -t sptree.tre Dsuite.vcf.gz SpeciesSet.txt
else
	Dsuite Dtrios -c -n gene_flow -t sptree.tre Dsuite.vcf.gz SpeciesSet.txt
fi
rm Dsuite.vcf.gz
echo -e "\nCreating Fbranch...."
if [[ $PBS_O_HOST == *".cz" ]]; then
	Dsuite/Build/Dsuite Fbranch sptree.tre SpeciesSet_gene_flow_tree.txt > SpeciesSet_gene_flow_Fbranch.txt
else
	Dsuite Fbranch sptree.tre SpeciesSet_gene_flow_tree.txt > SpeciesSet_gene_flow_Fbranch.txt
fi
if [[ $PBS_O_HOST == *".cz" ]]; then
	python3 Dsuite/utils/dtools.py -n gene_flow SpeciesSet_gene_flow_Fbranch.txt sptree.tre
else
	dtools.py -n gene_flow SpeciesSet_gene_flow_Fbranch.txt sptree.tre
fi
cairosvg gene_flow.svg -o gene_flow.pdf
echo -e "\nCreating Ruby files...\n"
cut -f 2 SpeciesSet.txt | uniq > plot_order.txt
ruby ./plot_d.rb SpeciesSet_gene_flow_BBAA.txt plot_order.txt 0.7 SpeciesSet_gene_flow_BBAA_D.svg
cairosvg SpeciesSet_gene_flow_BBAA_D.svg -o SpeciesSet_gene_flow_BBAA_D.pdf
ruby ./plot_f4ratio.rb SpeciesSet_gene_flow_BBAA.txt plot_order.txt 0.2 SpeciesSet_gene_flow_BBAA_f4ratio.svg
cairosvg SpeciesSet_gene_flow_BBAA_f4ratio.svg -o SpeciesSet_gene_flow_BBAA_f4ratio.pdf
rm plot_order.txt

#Copy results to home
if [[ $update =~ "yes" ]]; then
	cp *.{svg,png,pdf,txt,fasta,log,gz,tre} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Dsuite
else
	cp *.{svg,png,pdf,txt,fasta,log,gz,tre} $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Dsuite
fi

#Copy log to home
echo -e "\nEnd:" `date '+%A %d-%m-%Y %X'` >> ${logname}.log
if [[ $update =~ "yes" ]]; then
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/update/species_trees/Dsuite
else
	cp ${logname}.log $path/${treepath}${MISSINGPERCENT}_${SPECIESPRESENCE}/${tree}/species_trees/Dsuite
fi

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir08i
fi

echo -e "\nHybPhyloMaker8i finished...\n"
