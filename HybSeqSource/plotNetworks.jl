#!/usr/bin/env julia

#--------------------------------------------------------------------------------------------
# HybPhyloMaker: plotting phylogenetic networks using Julia
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0c
# Called from HybPhyloMaker8m2_PhyloNet_summary.sh
# Requires 'julia' and 'R'
# Install packages (PhyloNetworks, PhyloPlots, RCall) before running
# Tomas Fer, 2024
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------


# This script can be run as (example for outgroup defined in OUTGROUP):
# julia plotNetworks.jl "$OUTGROUP"
# without modifications plot phylogenetic networks in all files containing 'reti_' in their names
# one network per file, in rich newick format (e.g., results of PhyloNet or SNaQ)

outgroup = ARGS[1]
println("Outgroup: ", outgroup)

import Pkg
using PhyloNetworks
using PhyloPlots
using RCall

#define suffix for network plots
suffix = ".svg"
#define string for unrooted plots
unr = "_unrooted"
#define string for rooted plots
r = "_rooted"

#Loop over results
#get all the files containing 'reti_' in the current directory
files = filter(x -> occursin("reti_", x), readdir())
@info "Plotting unrooted networks"
for file in files
	#read the contents of the file into a variable with the same name as the file
	eval(Meta.parse("net = readTopology(\"$file\")")) #read the network from file
	imagefilename = "$file$unr$suffix"
	println(file)
	R"svg"(imagefilename, width=12, height=6)
	R"par"(mar=[0.1,0.1,0.1,0.1])
	plot(net, showgamma=true, tipcex=0.6, tipoffset = 0.2, style = :fulltree);
	R"dev.off()"
end

if outgroup != "test"
	#Plot rooted networks - this fails if rooting is incompatible with hybridization node!
	@info "Plotting rooted networks"
	for file in files
		#read the contents of the file into a variable with the same name as the file
		eval(Meta.parse("net = readTopology(\"$file\")")) #read the network from file
		imagefilename = "$file$r$suffix"
		println(file)
		#Reroot the network using outgroup
		try
			rootatnode!(net, outgroup)
			R"svg"(imagefilename, width=12, height=6)
			R"par"(mar=[0.1,0.1,0.1,0.1])
			plot(net, showgamma=true, tipcex=0.6, tipoffset = 0.2, style = :fulltree);
			R"dev.off()"
		catch e
			println("Caught an exception for $name: ", e)
		end
	end
else
	println("Outgroup was not set")
end

@info "Plotting done"
