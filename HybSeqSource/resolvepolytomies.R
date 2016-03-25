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
