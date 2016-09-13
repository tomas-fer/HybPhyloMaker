#!/usr/bin/python
#Combine bootstrap support values from two trees into one using p4 (Foster 2004)
#Modified from http://p4.nhm.ac.uk/tutorial/combine_supports.html
#modified for HybPhyloMaker by T. Fer, 2016
#-------------------------------------------------------------
#Two trees and alignment named 'combined.phylip' must be in the same directory
#Output tree 'combinedSupportsTree.tre' is in Newick format
#Usage combineboot.py masterTree secondaryTree

import sys

tree1 = sys.argv[1]
tree2 = sys.argv[2]
#alignment = sys.argv[3]

#print 'Tree1:', tree1
#print 'Tree2:', tree2
#print 'Alignment:', alignment

import p4
from p4 import *

#Disable checking for undetermined positions
var.doCheckForAllGapColumns=False
#Get a valid list of taxnames from an alignment (must be named 'concatenated.phylip'):
a = func.readAndPop('concatenated.phylip')
#print a.taxNames

#Read in and name the two trees, and make sure they both have the same taxNames:
tMaster = p4.func.readAndPop(tree1)
tSecondary = p4.func.readAndPop(tree2)

tMaster.taxNames = a.taxNames
tSecondary.taxNames = a.taxNames

# Split keys are numerical versions of the 'dot-star' split notation.
# The same split on the two trees would have the same split key.
tMaster.makeSplitKeys()
tSecondary.makeSplitKeys()

# Make a dictionary, so that we can fish out nodes in the secondary tree
# given a split key.  Split keys are found on node branches, here n.br.
myDict = {}
for n in tSecondary.iterInternalsNoRoot():
    myDict[n.br.splitKey] = n

for nM in tMaster.iterInternalsNoRoot():
    # Given a split key in the master tree, we can find the
    # corresponding node in the secondary tree, using the split key with
    # the dictionary.
    nS = myDict.get(nM.br.splitKey)
    # If there was none, then nS is None
    if nS:
        nM.name = '%s/%s' % (nM.name, nS.name)
    else:
        nM.name = '%s/-' % nM.name
    #print nM.name
tMaster.writeNewick('combinedSupportsTree.tre')

