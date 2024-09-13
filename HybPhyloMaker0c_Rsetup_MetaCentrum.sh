#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=4:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N Rpackages_setup
#PBS -m abe

# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *                Script 0c - Setup R-4.4 packages on Metacentrum               *
# *                                   v.1.8.0g                                   *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2024 *
# * tomas.fer@natur.cuni.cz                                                      *
# ********************************************************************************

#Download working versions of specific packages and install them to personal library
#ape 5.8 (necessary for phangorn 3.0.0)
#data.table 1.15.4
#phytools 2.1-1
#seqinr 4.2-36
#openxlsx 4.2.5.2
#gplots 3.1.3.1
#phangorn 3.0.0 (developmental version from GitHub)
#BiocManager 1.30.23 (for installation of Bioconductor packages like treeio)
#treeio 1.28.0 (from Bioconductor)
#all other package are required for smooth installation of these main packages

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
module add r/4.4.0-gcc-10.2.1-ssuwpvb

#Use this command before running R in HybPhyloMaker scripts
export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"

#Get R packages
echo -e "\nDownloading R packages (newest versions)...\n"
packages=( lattice nlme Rcpp digest ape data.table maps MASS clusterGeneration coda combinat iterators codetools foreach doParallel Matrix expm mnormt numDeriv optimParallel remotes BiocManager scatterplot3d fastmatch cli glue rlang lifecycle gtable isoband mgcv farver labeling colorspace munsell R6 RColorBrewer viridisLite scales fansi magrittr vctrs utf8 pillar tibble withr ggplot2 ggseqlogo igraph quadprog DEoptim phytools pixmap sp RcppArmadillo ade4 segmented seqinr stringi zip openxlsx tidyselect dplyr lazyeval fs fastmap cachem memoise yulab.utils purrr stringr tidyr tidytree gtools bitops caTools KernSmooth gplots )

for i in "${packages[@]}"; do
	version=$(wget -qO- https://cran.r-project.org/package=${i} | grep "tar.gz" | cut -d' ' -f6 | cut -d'_' -f2 | cut -d't' -f1)
	echo -e "${i}\t${version}"
	wget https://cran.r-project.org/src/contrib/${i}_${version}tar.gz
	echo -e "--------------------------------\n"
done

bcpackages=( treeio )
for i in "${bcpackages[@]}"; do
	version=$(wget -qO- https://bioconductor.org/packages/release/bioc/html/${i}.html | grep "tar.gz" | cut -d'_' -f3 | cut -d't' -f1)
	echo -e "${i}\t${version}"
	wget https://bioconductor.org/packages/release/bioc/src/contrib/${i}_${version}tar.gz
	echo -e "--------------------------------\n"
done

# #Direct links (might not work if newer version is released!)
# echo -e "\nDownloading R packages..."
# wget https://cran.r-project.org/src/contrib/lattice_0.22-6.tar.gz #required for ape
# wget https://cran.r-project.org/src/contrib/Archive/nlme/nlme_3.1-165.tar.gz #required for ape
# wget https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.12.tar.gz #required for ape
# wget https://cran.r-project.org/src/contrib/Archive/digest/digest_0.6.36.tar.gz #required for ape
# wget https://cran.r-project.org/src/contrib/ape_5.8.tar.gz
# wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.15.4.tar.gz
# wget https://cran.r-project.org/src/contrib/maps_3.4.2.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-60.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/clusterGeneration_1.3.8.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/coda_0.19-4.1.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/combinat_0.0-8.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/iterators_1.0.14.tar.gz #required for doParallel
# wget https://cran.r-project.org/src/contrib/codetools_0.2-20.tar.gz #required for foreach
# wget https://cran.r-project.org/src/contrib/foreach_1.5.2.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/doParallel_1.0.17.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/Matrix_1.7-0.tar.gz #required for expm
# wget https://cran.r-project.org/src/contrib/Archive/expm/expm_0.999-9.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/mnormt_2.1.1.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/numDeriv_2016.8-1.1.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/optimParallel_1.0-2.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/remotes_2.5.0.tar.gz
# wget https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.24.tar.gz
# wget https://klausvigo.r-universe.dev/bin/linux/noble/4.4/src/contrib/phangorn_3.0.0.0.tar.gz #this version not working on MetaCentrum, see below
# wget https://cran.r-project.org/src/contrib/scatterplot3d_0.3-44.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/fastmatch_1.1-4.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/cli_3.6.2.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/glue_1.7.0.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/rlang_1.1.4.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/lifecycle_1.0.4.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/gtable_0.3.5.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/isoband_0.2.7.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/mgcv_1.9-1.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/farver_2.1.2.tar.gz #required for scales
# wget https://cran.r-project.org/src/contrib/labeling_0.4.3.tar.gz #required for scales
# wget https://cran.r-project.org/src/contrib/colorspace_2.1-0.tar.gz #required for munsell
# wget https://cran.r-project.org/src/contrib/munsell_0.5.1.tar.gz #required for scales
# wget https://cran.r-project.org/src/contrib/R6_2.5.1.tar.gz #required for scales
# wget https://cran.r-project.org/src/contrib/RColorBrewer_1.1-3.tar.gz #required for scales
# wget https://cran.r-project.org/src/contrib/viridisLite_0.4.2.tar.gz #required for scales
# wget https://cran.r-project.org/src/contrib/scales_1.3.0.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/fansi_1.0.6.tar.gz #required for pillar
# wget https://cran.r-project.org/src/contrib/magrittr_2.0.3.tar.gz #required for pillar
# wget https://cran.r-project.org/src/contrib/vctrs_0.6.5.tar.gz #required for pillar
# wget https://cran.r-project.org/src/contrib/utf8_1.2.4.tar.gz  #required for pillar
# wget https://cran.r-project.org/src/contrib/pillar_1.9.0.tar.gz  #required for ggplot2
# wget https://cran.r-project.org/src/contrib/tibble_3.2.1.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/withr_3.0.0.tar.gz #required for ggplot2
# wget https://cran.r-project.org/src/contrib/ggplot2_3.5.1.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/ggseqlogo_0.2.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/igraph_2.0.3.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/quadprog_1.5-8.tar.gz #required for phytools
# wget https://cran.r-project.org/src/contrib/phytools_2.1-1.tar.gz
# wget https://cran.r-project.org/src/contrib/pixmap_0.4-13.tar.gz
# wget https://cran.r-project.org/src/contrib/sp_2.1-4.tar.gz
# wget https://cran.r-project.org/src/contrib/RcppArmadillo_0.12.8.4.0.tar.gz
# wget https://cran.r-project.org/src/contrib/ade4_1.7-22.tar.gz #required for seqinr
# wget https://cran.r-project.org/src/contrib/segmented_2.1-0.tar.gz #required for seqinr
# wget https://cran.r-project.org/src/contrib/seqinr_4.2-36.tar.gz
# wget https://cran.r-project.org/src/contrib/stringi_1.8.4.tar.gz #required for openxlsx
# wget https://cran.r-project.org/src/contrib/zip_2.3.1.tar.gz #required for openxlsx
# wget https://cran.r-project.org/src/contrib/openxlsx_4.2.5.2.tar.gz
# wget https://cran.r-project.org/src/contrib/tidyselect_1.2.1.tar.gz #required for dplyr
# wget https://cran.r-project.org/src/contrib/dplyr_1.1.4.tar.gz #required for treeio
# wget https://cran.r-project.org/src/contrib/lazyeval_0.2.2.tar.gz #required for tidytree
# wget https://cran.r-project.org/src/contrib/fs_1.6.4.tar.gz #required for yulab.utils
# wget https://cran.r-project.org/src/contrib/fastmap_1.2.0.tar.gz #required for cachem
# wget https://cran.r-project.org/src/contrib/cachem_1.1.0.tar.gz #required for memoise
# wget https://cran.r-project.org/src/contrib/memoise_2.0.1.tar.gz #required for yulab.utils
# wget https://cran.r-project.org/src/contrib/yulab.utils_0.1.7.tar.gz #required for treeio
# wget https://cran.r-project.org/src/contrib/purrr_1.0.2.tar.gz #required for tidyr
# wget https://cran.r-project.org/src/contrib/stringr_1.5.1.tar.gz #required for tidyr
# wget https://cran.r-project.org/src/contrib/tidyr_1.3.1.tar.gz #required for tidytree
# wget https://cran.r-project.org/src/contrib/tidytree_0.4.6.tar.gz #required for treeio
# wget https://cran.r-project.org/src/contrib/gtools_3.9.5.tar.gz #required for gplots
# wget https://cran.r-project.org/src/contrib/bitops_1.0-7.tar.gz #required for caTools
# wget https://cran.r-project.org/src/contrib/caTools_1.18.2.tar.gz #required for gplots
# wget https://cran.r-project.org/src/contrib/KernSmooth_2.23-24.tar.gz #required for gplots
# wget https://cran.r-project.org/src/contrib/gplots_3.1.3.1.tar.gz
# wget https://www.bioconductor.org/packages/release/bioc/src/contrib/treeio_1.28.0.tar.gz

#Install packages to personal (writable) library
mkdir -p /storage/$server/home/$LOGNAME/Rpackages44
packages=( lattice nlme Rcpp digest ape data.table maps MASS clusterGeneration coda combinat iterators codetools foreach doParallel Matrix expm mnormt numDeriv optimParallel remotes BiocManager )
for package in "${packages[@]}"; do
	echo -e "\nInstalling $package..."
	R CMD INSTALL ${package}_* -l /storage/$server/home/$LOGNAME/Rpackages44
done

#Install developmental version of phangorn 3.0.0 from GitHub
R -q -e "library(BiocManager);BiocManager::install('Biostrings')"
R -q -e "library(remotes);remotes::install_github('KlausVigo/phangorn')"

packages2=( scatterplot3d fastmatch cli glue rlang lifecycle gtable isoband mgcv farver labeling colorspace munsell R6 RColorBrewer viridisLite scales fansi magrittr vctrs utf8 pillar tibble withr ggplot2 ggseqlogo igraph quadprog DEoptim phytools pixmap sp RcppArmadillo ade4 segmented seqinr stringi zip openxlsx tidyselect dplyr lazyeval fs fastmap cachem memoise yulab.utils purrr stringr tidyr tidytree treeio gtools bitops caTools KernSmooth gplots )

for package in "${packages2[@]}"; do
	echo -e "\nInstalling $package..."
	R CMD INSTALL ${package}* -l /storage/$server/home/$LOGNAME/Rpackages44
done

#Check R packages
echo -e "\n**************************************************************"
echo -e "Checking R packages"
for Rpackage in ape seqinr data.table openxlsx phytools phangorn treeio gplots; do
	R -q -e "aa <- file('Rtest', open='wt'); sink(aa, type='message'); require($Rpackage); sink(type='message'); close(aa)" > /dev/null
	if grep -Fq "no package called" Rtest; then
		echo -e "R package $Rpackage...not found"
	elif grep -Fq "Error" Rtest; then
		echo -e "R package $Rpackage...unable to load"
	else
		ver=$(R -q -e "packageVersion('$Rpackage')" 2> /dev/null | grep "[1]" | cut -d' ' -f2 | sed 's/^.//' | sed 's/.$//')
		echo -e "R package $Rpackage...OK...version: ${ver}"
	fi
done
rm Rtest
echo -e "**************************************************************\n"

#Delete scratch
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
fi

echo -e "\nScript HybPhyloMaker0c finished...\n"
