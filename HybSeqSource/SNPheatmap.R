#--------------------------------------------------------------------------------------------
# HybPhyloMaker: similarity heatmap based on filtered SNPs from concatenated alignment
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0a
# Tomas Fer, 2025
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------

library(gplots) #for heatmap

#Read arguments from command line
args <- commandArgs()
calctype <- as.numeric(args[5])

#Read data matrix
a <- read.table("SNPmatrixTransposed.txt", header=T)

#Calculate pair-wise distance matrix
if (calctype == 0) {
  #Calculate 'simple matching' coefficient using PERMANOVA package
  library (PERMANOVA)
  mat <- a[, -1]
  names <- a[, 1]
  distbin <- DistBinary(mat, coefficient = "Simple_Matching", transformation = "1-S")
  distmatrix <- distbin$D
} else {
  #Calculate similarity coefficient as a percentage of matches without NA fields
  library(dplyr)
  #Define simple matching coefficient function
  simple_matching_coefficient <- function(s1, s2) {
    #Remove NA values
    notna <- !is.na(s1) & !is.na(s2)
    s1 <- s1[notna]
    s2 <- s2[notna]
    #Calculate matching coefficient
    matches <- sum(s1 == s2)
    total <- length(s1)
    return(matches/total)
  }
  nsamples <- nrow(a)
  names <- a$POS
  #Prepare empty matrix
  similmatrix <- matrix(NA, nrow = nsamples, ncol = nsamples)
  rownames(similmatrix) <- names
  colnames(similmatrix) <- names
  #Nested loop over all samples versus all samples
  for (i in 1:nsamples) {
    for (j in i:nsamples) {
      sample1 <- a[i, -1] %>% unlist()
      sample2 <- a[j, -1] %>% unlist()
      similmatrix[i, j] <- simple_matching_coefficient(sample1, sample2)
      similmatrix[j, i] <- similmatrix[i, j] #symmetric matrix
    }
  }
  #Make distance coefficient
  distmatrix <- 1-similmatrix
}

#Plot&save heatmap
#colors customized from 'scico' lajolla palette (8 colours)
#set color palette (first color(black) is given 4x to increase its length in gradient
colpal=c("#191900", "#191900", "#191900", "#191900", "#191900", "#3C2614", "#773831", "#C04D49", "#E0764F", "#EAA252", "#F5D268", "#F8DE7A", "#FFFECB")
pdf(file="SNPheatmap.pdf", width=12, height=12)
hm <- heatmap.2(distmatrix, density.info=c("none"), keysize=0.5, key.title="Dissimilarity", key.xlab=NA, lwid=c(1, 7), lhei=c(1, 13), trace="none", Rowv=T, Colv=T, cexRow=0.8, cexCol=0.8, labRow=names, labCol=names, margins=c(14, 14), col=colorRampPalette(colpal)(170))
dev.off()

#Save distance matrix to file
distmat <- hm$carpet
colnames(distmat) <- names[hm$colInd]
rownames(distmat) <- names[hm$rowInd]
write.table(distmat, 'distmat.txt', col.names=NA, sep="\t")
