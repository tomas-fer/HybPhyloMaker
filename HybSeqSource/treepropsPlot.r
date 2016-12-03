#MAKE PNG PICTURE WITH HISTOGRAM AND BOXPLOT FOR TREE PROPERTIES
#based on http://rgraphgallery.blogspot.com/2013/04/rg-plotting-boxplot-and-histogram.html
#modified for HybPhyloMaker by T. Fer, 2016
#-------------------------------------------------------------

#Read table from 'tree_stats_table.csv' 
x <- read.csv("tree_stats_table.csv", header=TRUE, stringsAsFactors=F)

#Make PNG picture of histogram from values in a table (omitting first column with names)
for (i in 2:length(colnames(x))) {
  png(file=paste(colnames(x[i]), "_histogram.png", sep = ""), width=600, height=600)
  nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))
  par(mar=c(3.1, 3.1, 1.1, 2.1))
  boxplot(x[,i], horizontal=TRUE,  outline=TRUE, ylim=c(min(x[,i]), max(x[,i])), frame=F, axes = FALSE, col = "grey")
  hist(x[,i],xlim=c(min(x[,i]), max(x[,i])), col = "lightblue", main = colnames(x[i]))
  dev.off()
}
