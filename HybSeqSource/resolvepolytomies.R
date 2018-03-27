#--------------------------------------------------------------------------------------------
# HybPhyloMaker: resolve polytomies using multi2di 
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.6.0
# Taken from K. M. Everson's blog (https://www.kmeverson.org/blog/force-bifurcating-trees-in-r)
# Modified for HybPhyloMaker
# Tomas Fer, 2018
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------

library(ape)

args <- commandArgs()
name <- args[5]

myTrees<-read.tree(name)

make.bifurcating<-function(phy){
  for(i in 1:length(phy)){
    phy[[i]]<-multi2di(phy[[i]])
  }
  return(phy)
}

myTrees.bifurcating<-make.bifurcating(myTrees)

#is.binary.tree(myTrees.bifurcating[[61]])

write.tree(myTrees.bifurcating, file=paste0(name,".bifurcating.newick"))
