#--------------------------------------------------------------------------------------------
# HybPhyloMaker: similarity heatmap based on filtered SNPs from concatenated alignment
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0
# Tomas Fer, 2025
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------

library (PERMANOVA) #for distance matrix
library(gplots) #for heatmap

#Read and parse data matrix
a <- read.table("SNPmatrixTransposed.txt", header=T)
mat <- a[,-1]
names <- a[,1]

#Calculate pair-wise distance matrix using 'simple matching' coefficient
distbin <- DistBinary(mat, coefficient = "Simple_Matching")

#Plot&save heatmap
#colors customized from 'scico' lajolla palette (8 colours)
#set color palette (first color(black) is given 4x to increase its length in gradient
colpal=c("#191900","#191900","#191900","#191900","#191900","#3C2614","#773831","#C04D49","#E0764F","#EAA252","#F5D268","#F8DE7A","#FFFECB")
pdf(file="SNPheatmap.pdf", width=12, height=12)
hm <- heatmap.2(distbin$D,density.info=c("none"),keysize=0.5,key.title="Dissimilarity",key.xlab=NA,lwid=c(1,7),lhei=c(1,13),trace="none",Rowv=T,Colv=T,cexRow=0.8,cexCol=0.8,labRow=names,labCol=names,margins=c(14, 14),col=colorRampPalette(colpal)(170))
dev.off()

#Save distance matrix to file
distmat <- hm$carpet
colnames(distmat) <- names[hm$colInd]
rownames(distmat) <- names[hm$rowInd]
write.table(distmat, 'distmat.txt', col.names=NA, sep = "\t")
