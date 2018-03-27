#--------------------------------------------------------------------------------------------
# HybPhyloMaker: calculating tree and alignment properties
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.6.0
# Taken from M. Borowiec GitHub (https://github.com/marekborowiec/good_genes)
# Modified for HybPhyloMaker
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
#all alignments must be in dir 'alignments'
aln_dir <- file.path("./alignments/")

#names of the trees must be 'something_LOCUSNR_somethingelse.tre'
trees_files <- dir(path=trees_dir, pattern="*tre")
#alignment files must be 'something.fas'
aln_files <- dir(path=aln_dir, pattern="*fas")

### AVERAGE BOOTSTRAP SUPPORT ###
# This code should work with any Newick tree files
# that have some measure of support at nodes, 
# including PP from Bayesian analysis.

# define a function to calculate average support
Avg_support <- function(file) {
  
  # get Assembly number from filename (second part of file name separated by '_')
  locus_no <- strsplit(file, "_")
  locus_no <- locus_no[[1]][2]
  # read in tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # store support values in a vector
  support <- c(as.numeric(tree$node.label))
  # calculate average support
  avg_supp <- mean(support, na.rm=T)
  return(c(locus_no,avg_supp))
  
}

# loop over all files
average_bootstrap <- lapply(trees_files, Avg_support)
average_bootstrap <- data.frame(matrix(unlist(average_bootstrap), nrow=(length(average_bootstrap)), byrow=T))
colnames(average_bootstrap) <- c("Locus", "Average_bootstrap")


### AVERAGE BRANCH LENGTHS ###
# This takes a Newick tree with branch lengths
# and returns the number of tips for each tree,
# and calculates average branch length.

Br_length.trees <- function(file) {
  
  # reads the phylogenetic tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # gets number of tips
  no_tips <- length(tree$tip.label)
  # calculate avg branch length
  avg_br_length <- mean(tree$edge.length)
  # get Assembly number from filename (second part of file name separated by '_')
  locus_no <- strsplit(file, "_")
  locus_no <- locus_no[[1]][2]
  return(c(locus_no,avg_br_length))
  
}

# loop over all files
br_lengths <- lapply(trees_files, Br_length.trees)
br_lengths <- data.frame(matrix(unlist(br_lengths), nrow=(length(br_lengths)), byrow=T))
colnames(br_lengths) <- c("Locus", "Average_branch_length")


### CLOCKLIKENESS ###
# define a function to calculate coefficient of variation
# of root to tip distances for the root minimizing the coefficient

Clocklikeness <- function(file) {
  
  # get Assembly number from filename (second part of file name separated by '_')
  locus_no <- strsplit(file, "_")
  locus_no <- locus_no[[1]][2]
  # read in tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # record coefficient for all possible outgroups
  CV <- c()
  taxa <- tree$tip.label
  for (taxon in taxa) {
    # root tree
    rooted_tr <- root(phy=tree,outgroup=taxon,resolve.root=T)
    # get matrix diagonal of phylogenetic variance-covariance matrix
    # these are your distances from root
    root_dist <- diag(vcv.phylo(rooted_tr))
    std_dev_root_dist <- sd(root_dist)
    mean_root_dist <- mean(root_dist)
    CV <- c(CV, (std_dev_root_dist/mean_root_dist)*100)
  }
  # get lowest CV
  minCV <- min(CV)
  print(minCV)
  return(c(locus_no,minCV))
  
}

# loop over all files
cv_clocklikeness <- lapply(trees_files, Clocklikeness)
cv_clocklikeness <- data.frame(matrix(unlist(cv_clocklikeness), nrow=(length(cv_clocklikeness)), byrow=T))
colnames(cv_clocklikeness) <- c("Locus", "Clocklikeness")


### SATURATION ###
# The code below takes as input FASTA alignments and tree files with branch lengths
# and calculates regression slope and r-squared for each locus
# also saving regression plots for each locus 
# in a newly created 'Saturation_plots' folder

# define function to perform regression, 
# calculate R-squared, and save saturation plots
# needs fasta alignments and 
# corresponding raxml's tree files
Saturation <- function(seq, tree, locus_no) {
  
  # read in alignment
  alignment <- read.alignment(file=seq, format="fasta")
  # read in tree
  tree <- read.tree(tree)
  # get locus number from filename
  # this has to be different regex from other functions used here
  # because file names here include directory
  locus_no <- strsplit(locus_no, "_")
  locus_no <- locus_no[[1]][2]
  # matrix with pairwise identity
  mat <- dist.alignment(alignment, matrix="identity")
  # matrix with uncorrected p-distances
  p_mat <- mat*mat
  # mean p-distance
  avg_p_mat <- mean(p_mat, na.rm=TRUE)
  # make matrix of pairwise distances in branch lengths from the tree
  cophentr <- cophenetic(tree)  
  # store as matrix objects
  mat_mat <- as.matrix(mat)
  mat_p_mat <- as.matrix(p_mat)
  # order p-distance matrix by names
  mat_p_mat <- mat_p_mat[order(row.names(mat_p_mat)),order(row.names(mat_p_mat))]
  mat_co <- as.matrix(cophentr)
  # order pairwise distances matrix by names
  mat_co <- mat_co[order(row.names(mat_co)),order(row.names(mat_co))]
  # get lower triangulars of both matrices
  branch_dist <- mat_co[lower.tri(mat_co)]
  p_dist <- mat_p_mat[lower.tri(mat_p_mat)]
  # perform simple linear regression
  regress <- lm(p_dist ~ branch_dist)
  # get slope
  slope <- coef(regress)[2]
  # get r-squared
  Rsquared <- summary(regress)$r.squared
  
  # plot branch length pairwise distances on x
  # and uncorrected p-distances on y
  
  # open png file
  #png(file=paste(sat_dir, locus_no, "-saturation.png", sep=""), width=600, height=600)
  #plot(branch_dist, p_dist)
  # add simple linear regression line
  #abline(lm(p_dist ~ branch_dist), col="red")
  # give title as locus number and subtitle as tree length
  #title(main=locus_no,sub=paste("Slope: ", round(slope, digits=3), " R-squared: ", round(Rsquared, digits=3), sep=""), cex.sub=1.25)
  # close png file
  #dev.off()
  
  return(list(locus_no, avg_p_mat, slope, Rsquared))
  
}

# get current working directory
#work_dir <- getwd()
# create a folder for saturation plots
#dir.create("./saturation_plots")
#sat_dir <- file.path(work_dir, "saturation_plots/")

# create a table with file names and column for filtering locus name
files_table <- as.data.frame(cbind(paste(aln_dir, aln_files, sep=""), paste(trees_dir, trees_files, sep=""), trees_files))
saturation_table <- t(mapply(Saturation, as.matrix(files_table$V1), as.matrix(files_table$V2), as.matrix(files_table$trees_files)))
saturation_table <- as.data.frame(saturation_table)
colnames(saturation_table) <- c("Locus","Avg_p_dist","Slope","R_squared")
row.names(saturation_table) <- NULL


### PLOTTING GENE TREES ### 
# This will save png files of 600 x 600 pixel plots
# of all unrooted phylogenetic trees with tip labels
# and support values

# Plot_trees <- function(file) {
#   
#   # reads the phylogenetic tree
#   tree <- read.tree(paste(trees_dir,file, sep=""))
#   # extracts plot name (locus) from file name 
#   plot_name <- sub(tree_regex, "\\1", perl=TRUE, x=file)
#   plot_name <- sub("Assembly_","", plot_name, fixed=TRUE)
#   # open png file
#   png(file=paste(tree_plots_dir, plot_name, "-tree.png", sep=""), width=600, height=1200)
#   plot.phylo(tree, show.node.label=T)
#   # give title as locus number and subtitle as tree length
#   title(main=plot_name)
#   # close png file
#   dev.off()
#   
# }

# # create directory for tree plots
# dir.create("./tree_plots")
# tree_plots_dir <- file.path("tree_plots/")
# 
# #loop over all files
# lapply(trees_files, Plot_trees)

# putting together all the data
dfs <- list(average_bootstrap, br_lengths, cv_clocklikeness, saturation_table)

Multmerge <- function(dfs){
  datalist <- lapply(dfs, function(x){data.frame(x)})
  Reduce(function(x,y) {merge(x,y)}, datalist)
}

all_loci_stats <- Multmerge(dfs)
all_loci_stats <- data.frame(lapply(all_loci_stats, as.character), stringsAsFactors=FALSE)

write.csv(all_loci_stats, file="tree_stats_table.csv", quote=FALSE, row.names=FALSE)
