#--------------------------------------------------------------------------------------------
# HybPhyloMaker: plot heatmap of missing data per sample and locus (for selected loci)
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0
# Tomas Fer, 2021
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------

#Processing table with missing data
#Creates a PDF image showing heatmap using heatmap.2 from 'gplots'
#(scaling of margins and label sizes might not work for all cases)

library(gplots)

#Read the data and remove unwanted columns and rows
data <- read.csv ("MissingDataOverview.txt", sep=" ")
data<-subset(data, species!="average_missing" & species!="total_genes")
maxspecies <- max(nchar(names(data)))
maxlocus <- max(nchar(as.vector(data$species)))
nrspecies <- length(names(data))
if(nrspecies < 24){nrspecies = 24}
nrloci <- length(data$species)
row.names(data) <- data$species
data$species <- NULL
data$averageMissing <- NULL
data$percPresentSpecies <- NULL
data_matrix <- data.matrix(data)

#Plot the heatmap (modify margins and cexCol/cexRow according to your needs)
pdf("heatmap.pdf")
heatmap.2(data_matrix, col = gray.colors(256), margins=c(3+(maxspecies*0.3),1+(maxlocus*0.3)), cexCol=0.2*((1/nrspecies)*100), cexRow=0.5, tracecol=NA, density.info="none", key.title=NA, key=F)
dev.off()
