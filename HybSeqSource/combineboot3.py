#!/usr/bin/python
#-----------------------------------------------------------------------------------------------
# HybPhyloMaker: combine bootstrap support values from two trees into one using p4 (Foster 2004)
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0
# Taken from http://p4.nhm.ac.uk/tutorial/combine_supports.html
# Modified for HybPhyloMaker
# Tomas Fer, 2025
# tomas.fer@natur.cuni.cz
#-----------------------------------------------------------------------------------------------

#Two trees and alignment named 'concatenated.phylip' must be in the same directory
#Output tree 'combinedSupportsTree.nex' is in NEXUS format
#Usage: python3 ./combineboot3.py tree1 tree2


import sys

print ('Number of arguments:', len(sys.argv))
tree1 = sys.argv[1]
tree2 = sys.argv[2]
#alignment = sys.argv[3]

print ('Tree1:', tree1)
print ('Tree2:', tree2)
#print 'Alignment:', alignment

import p4
from p4 import *

#we need a valid list of taxnames, which we might get from an alignment:
var.doCheckForAllGapColumns=False
# read(alignment)
# a = var.alignments[0]
a = func.readAndPop('concatenated.phylip')

print ()
print ('Taxnames')
print (a.taxNames)
print ()

#We read in and name our two trees, and make sure they both have the same taxNames:
tMaster = p4.func.readAndPop(tree1)
tSecondary = p4.func.readAndPop(tree2)

print ('Master tree')
tMaster.draw()
print ()
print ('Secondary tree')
tSecondary.draw()

tMaster.taxNames = a.taxNames
tSecondary.taxNames = a.taxNames

# Split keys are numerical versions of the 'dot-star' split notation.
# The same split on the two trees would have the same split key.
tMaster.makeSplitKeys()
tSecondary.makeSplitKeys()

# Make a dictionary, so that we can fish out nodes in the secondary tree
# given a split key.  Split keys are found on node branches, here
# n.br.
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
tMaster.writeNexus('combinedSupportsTree.nex')

