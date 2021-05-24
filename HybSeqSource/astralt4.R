#--------------------------------------------------------------------------------------------
# HybPhyloMaker: draw phylogenetic tree with pie charts based on 'ASTRAL -t 4'
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0
# Tomas Fer, 2021
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------

#Processing scored tree from ASTRAL (with the option -t 4)
#This creates a PDF image of the rooted ASTRAL species tree with pie charts on branches
#See ASTRAL manual for explanation

require(treeio) #(implements 'read.astral')
require(ape)

#Import the tree generated with ASTRAL -t 4, i.e. with three posterior probabilities per branch
x <- read.astral("tree.tre")
#Save posterior probabilities 
pp1 <- as.numeric(x[['pp1']])
pp2 <- as.numeric(x[['pp2']])
pp3 <- as.numeric(x[['pp3']])
pc <- cbind( pp1, pp2, pp3)

tree <- x@phylo

#Plot/save the tree (modify 'adj' and 'cex' values according to your needs)
pdf("tree.pdf")
par(lty = "blank") #no lines in pie charts
scaling = 0.27*((1/length(tree$tip.label))*100) #this should automatically scale tip label according to the number of tips
if(scaling > 1){scaling = 1} # set scaling to '1' if too few tips
plot(tree, cex = scaling, label.offset = 1)
nodelabels(pie = pc, adj = c(-0.2, 0.5), piecol = c("blue","red","green"), cex = 0.4)
dev.off()
