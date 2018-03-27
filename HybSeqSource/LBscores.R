#--------------------------------------------------------------------------------------------
# HybPhyloMaker: calculating of the long-branch scores
# https://github.com/tomas-fer/HybPhyloMaker
# Taken from M. Borowiec GitHub (https://github.com/marekborowiec/metazoan_phylogenomics)
# Modified for HybPhyloMaker
# v.1.6.0
# Tomas Fer, 2018
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------

# disable scientific notation
options(scipen=999)

# install required libraries
# uncomment if you don't have these installed
#install.packages("ape")
#install.packages("seqinr")
#install.packages("data.table")

# load needed libraries
library("ape")
library("seqinr")
library("data.table")

### SET WORKING DIRECTORY ###
#all trees must be in dir 'trees'
trees_dir <- file.path("./trees/")

#names of the trees must be 'something_LOCUSNR_somethingelse.tre'
trees_files <- dir(path=trees_dir, pattern="*tre")

### LONG-BRANCH SCORES ###
# This code calculates long-branch scores 
# as defined by Struck 2014 Evol. Bioinform. 10:51-67
# for each taxon in each locus, as well as standard deviation
# of these scores in each locus

LB_score <- function(file) {
  
  # read the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # make matrix of pairwise distances in branch lengths from the tree
  cophentr <- cophenetic(tree)
  # get Assembly number from filename (second part of file name separated by '_')
  locus_no <- strsplit(file, "_")
  locus_no <- locus_no[[1]][2]
  # create empty data frame
  LB_table <- list()
  
  # loop over all taxa in the matrix
  Scoring <- function(taxon) {
    
    # calculate mean of patristic distances for the taxon
    tax_mean <- mean(cophentr[,taxon])
    # calculate Struck's LB score
    LB_score <- ( tax_mean / mean(cophentr) - 1 ) * 100
    # define a row for the table
    return <- c(taxon, tax_mean, LB_score, locus_no)
    
  }
  
  LB_table <- lapply(row.names(cophentr), Scoring)
  LB_table <- data.frame(matrix(unlist(LB_table), nrow=(length(LB_table)), byrow=T))
  colnames(LB_table) <- c("Taxon", "Tax_mean", "LB_score", "Locus")
  
  # calculate standard deviation of LB
  LB_SD <- sd(LB_table$LB_score, na.rm=T)
  return(data.frame(LB_SD, LB_table))
  
}

# loop over all files
LB_scores <- lapply(trees_files, LB_score)
# compress the list to a data frame
LB_scores <- rbindlist(LB_scores)

write.csv(LB_scores, file="LBscores.csv", quote=FALSE, row.names=FALSE)