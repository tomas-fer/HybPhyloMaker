#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2h
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l mem=1gb
#PBS -N HybPipe0_preparereference
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N HybPipe0_preparereference
#$ -o HybPipe0_preparereference.log

# ********************************************************************************
# *       HybPipe - Pipeline for Hyb-Seq data processing and tree building       *
# *                     Script 00b - Prepare pseudoreference                     *
# *                                   v.1.0.0                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2016 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

# Make pseudo-reference from multifasta file (e.g. Hyb-Seq probes), i.e. make single fasta (one line) with all original sequences separated by 400 Ns
# 400 Ns are also added to the beginning and the end of the reference
# Pseudo-reference (*_with400Ns_beginend.fas) is saved to HybSeqSource directory

if [[ $PBS_O_HOST == *".cz" ]]; then
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
elif [[ $HOSTNAME == *local* ]]; then
	echo "Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir00_ref
	cd workdir00_ref
else
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir workdir00_ref
	cd workdir00_ref
fi

#Copy probe sequence file
cp $source/$probes .

# Cut filename before first '.', i.e. remove suffix - does not work if there are other dots in reference file name
name=`ls $probes | cut -d'.' -f 1`
# Print headers in fasta file
echo ">${name}_with400Ns_beginend" > ${name}_with400Ns_beginend.fas
# 0. print N 400x to variable $a
a=$(printf "%0.sN" {1..400})
# 1. awk command to remove EOLs from lines not beginning with '>', i.e. all lines containing sequence for particular record are merged, i.e. each record on two lines only
# 2. sed command to delete all lines starting with '>', i.e. only sequences remain
# 3. tr command to replace all EOLs ("\n") by Q
# 4. sed command to replace all Qs by a sequence of 400 Ns (saved in variable $a)
# 5. awk command to print 400 Ns to the beginning and the end of the reference
cat $probes | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed '/>/d' | tr "\n" "Q" | sed "s/Q/$a/g" | awk -v val=$a '{ print val $0 val }' >> ${name}_with400Ns_beginend.fas

#Copy pseudoreference to home
cp ${name}_with400Ns_beginend.fas $source

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	rm -rf $SCRATCHDIR/*
else
	cd ..
	rm -r workdir00_ref
fi
