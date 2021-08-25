#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N Rpackages_setup
#PBS -m abe

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                 Script 0c - Setup R packages on Metacentrum                  *
# *                                   v.1.8.0a                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2021 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Download working versions of specific packages and install them to personal library
#ape 5.0 (IMPORTANT, v 3.4 is not compatible with HybPhyloMaker)
#data.table 1.9.6
#phytools 0.5-20
#seqinr 3.1-3
#openxlsx 4.0.0
#phangorn 2.5.5 (this works also with R 3.5 and lower)
#BiocManager 1.30.14 (for installation of Bioconductor packages like treeio)
#treeio 1.4.3 (from Bioconductor, this works with R 3.4)
#gplots 3.0.1.2 (this works with R 3.4)

if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nScript HybPhyloMaker0c is running on Metacentrum..."
else
	echo -e "\nYou are not on Metacentrum. Install the R packages using 'install_software.sh' script. Exiting...\n"
	exit 3
fi

#Move to scratch
cd $SCRATCHDIR
#Copy file with settings from home and set variables from settings.cfg
cp $PBS_O_WORKDIR/settings.cfg .
. settings.cfg

#Add necessary modules
module add R-3.4.3-gcc

#Use this command before running R in HybPhyloMaker scripts
export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages"

#Get ape 3.5 and other R packages
echo -e "\nDownloading R packages..."
wget https://cran.r-project.org/src/contrib/Archive/ape/ape_5.0.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.9.6.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/phytools/phytools_0.5-20.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/seqinr/seqinr_3.1-3.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/openxlsx/openxlsx_4.0.0.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/quadprog/quadprog_1.5-7.tar.gz #required by phangorn
wget https://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.2.5.tar.gz #required by phangorn
wget https://cran.r-project.org/src/contrib/fastmatch_1.1-0.tar.gz #required by phangorn
wget https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.14.tar.gz #necessary for treeio on Bioconductor
wget https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.7.tar.gz #required by treeio
wget https://cran.r-project.org/src/contrib/Archive/tidytree/tidytree_0.3.3.tar.gz #required by treeio
wget https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz #required by treeio
wget https://cran.r-project.org/src/contrib/Archive/bitops/bitops_1.0-6.tar.gz #required by gplots
wget https://cran.r-project.org/src/contrib/Archive/gdata/gdata_2.17.0.tar.gz #required by gplots
wget https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.4.tar.gz #required by gplots

#Install packages to personal (writable) library
mkdir -p /storage/$server/home/$LOGNAME/Rpackages
for package in $(ls *.tar.gz); do
	echo -e "\nInstalling $package..."
	R CMD INSTALL $package -l /storage/$server/home/$LOGNAME/Rpackages
done

#instal 'phangorn'
wget https://cran.r-project.org/src/contrib/Archive/phangorn/phangorn_2.5.5.tar.gz
R CMD INSTALL phangorn_2.5.5.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #this only works once the previous packages are installed

#install 'ellipsis'
wget https://cran.r-project.org/src/contrib/Archive/ellipsis/ellipsis_0.3.1.tar.gz
R CMD INSTALL ellipsis_0.3.1.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #required by dplyr

#install 'Rcpp'
wget https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.1.tar.gz
R CMD INSTALL Rcpp_1.0.1.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #required by dplyr

#install 'glue'
wget https://cran.r-project.org/src/contrib/Archive/glue/glue_1.4.1.tar.gz
R CMD INSTALL glue_1.4.1.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #required by dplyr

#install 'pillar'
wget https://cran.r-project.org/src/contrib/Archive/pillar/pillar_1.3.1.tar.gz
R CMD INSTALL pillar_1.3.1.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #required by tibble

#install 'tibble'
wget https://cran.r-project.org/src/contrib/Archive/tibble/tibble_2.1.3.tar.gz
R CMD INSTALL tibble_2.1.3.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #required by dplyr

#install 'tidyselect'
wget https://cran.r-project.org/src/contrib/Archive/tidyselect/tidyselect_0.2.5.tar.gz
R CMD INSTALL tidyselect_0.2.5.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #required by dplyr

#install 'dplyr'
wget https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.8.5.tar.gz #required by rlang
R CMD INSTALL dplyr_0.8.5.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages

#instal 'treeio'
wget https://www.bioconductor.org/packages/3.7/bioc/src/contrib/treeio_1.4.3.tar.gz
R CMD INSTALL treeio_1.4.3.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #this only works once the previous packages are installed

#instal 'gplots'
wget https://cran.r-project.org/src/contrib/Archive/gplots/gplots_3.0.1.2.tar.gz
R CMD INSTALL gplots_3.0.1.2.tar.gz -l /storage/$server/home/$LOGNAME/Rpackages #this only works once the previous packages are installed

#Delete scratch
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
fi

echo -e "\nScript HybPhyloMaker0c finished...\n"
