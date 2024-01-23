#!/usr/bin/env julia

#--------------------------------------------------------------------------------------------
# HybPhyloMaker: running SNaQ analysis on a set of gene trees
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0c
# Called from HybPhyloMaker8l_SNaQ.sh
# Requires 'julia' and 'R'
# Install packages (PhyloNetworks, PhyloPlots, RCall, DataFrames, CSV, Distributed) before running
# Cite 'PhyloNetworks' (https://github.com/crsl4/PhyloNetworks.jl) when using this script!
# Implemented by
# Roman Ufimov & Tomas Fer, 2023
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------


# This script can be run as (example for outgroup defined in OUTGROUP, hmin=0, hmax=5 and 10 runs):
# julia runSNaQ.jl "$OUTGROUP" 0 5 10

length(ARGS) > 2 ||
    error("Julia: Need 3 or 4 arguments: outgroup, min nr reticulations, max nr reticulations and nr runs (optional, 10 by default)")
outgroup = ARGS[1]
h_start = parse(Int, ARGS[2])
h_max = parse(Int, ARGS[3])
nruns = 10
println("Outgroup: ", outgroup, "\nhmin: ", h_start, "\nhmax: ", h_max)
if length(ARGS) > 3
    nruns = parse(Int, ARGS[4])
end
println("nruns: ", nruns)

import Pkg
using PhyloNetworks
using PhyloPlots
using RCall
using DataFrames
using CSV
using Distributed
addprocs(nruns)
@everywhere using PhyloNetworks

@info "Creating CF file"
genetrees = readMultiTopology("trees.newick");
q,t = countquartetsintrees(genetrees);
df = writeTableCF(q,t)
using CSV
CSV.write("tableCF.csv", df);
raxmlCF = readTableCF(df);

h=h_start
h_before = h - 1
if h_before < 0
    @info "Reading species tre"
    tre = readTopology("sptree.tre")
    outputfile = string("net", h, "_", nruns, "runs") # example: "net2_10runs"
    seed = 1234 + h # change as desired! Best to have it different for different h
    @info "Will run SNaQ with h=$h, # of runs=$nruns, seed=$seed, output will go to: $outputfile"
    @info "Running SNaQ..."
    net0 = snaq!(tre, raxmlCF, hmax=h, filename=outputfile, seed=seed, runs=nruns)
    net = net0
else
    net = readSnaqNetwork(string("net", h_before, "_", nruns, "runs.out"))
end

while h <= h_max
    outputfile = string("net", h, "_", nruns, "runs") # example: "net2_10runs"
    seed = 1234 + h # change as desired! Best to have it different for different h
    @info "Will run SNaQ with h=$h, # of runs=$nruns, seed=$seed, output will go to: $outputfile"
    @info "Running SNaQ"
    current_net = snaq!(net, raxmlCF, hmax=h, filename=outputfile, seed=seed, runs=nruns)
    global net = current_net
    global h = h + 1 
end

@info "Making summary plots"
#define suffix for network plots
suffix = ".svg"
#define string for unrooted plots
unr = "_unrooted"
#define string for rooted plots
r = "_rooted"
#define empty vector for network likelihood scores
scores = Vector{Float64}()

#Loop over results
#get all the files with the .out suffix in the current directory
files = filter(x -> occursin(".out", x), readdir())
@info "Plotting unrooted networks"
for file in files
	#extract the filename without the extension
	name = splitext(file)[1]
	#read the contents of the file into a variable with the same name as the file
	eval(Meta.parse("const $name = readSnaqNetwork(\"$file\")")) #read the whole output
	#put network to variable 'net'
	eval(Meta.parse("net = $name"))
	#put network likelihood to variable 'netlik'
	eval(Meta.parse("netlik = $name.loglik"))
	#add likelihood to a vector 'scores'
	push!(scores, netlik)
	imagefilename = "$name$unr$suffix"
	println(name)
	R"svg"(imagefilename, width=12, height=6)
	R"par"(mar=[0.1,0.1,0.1,0.1])
	plot(net, showgamma=true, tipcex=0.6, tipoffset = 0.2, style = :fulltree);
	R"dev.off()"
end

#Plot scores
@info "Plotting network scores"
hm = collect(0:length(files)-1)
R"plot"(hm, scores, type="b", ylab="Network score", xlab="Hmax", col="blue");
R"dev.off()"
#Save network scores to a file
y = DataFrame(Hmax = hm, NetScore = scores)
CSV.write("scores.txt", y; delim='\t')

#Plot rooted networks - this fails if rooting is incompatible with hybridization node!
@info "Plotting rooted networks"
for file in files
	#extract the filename without the extension
	name = splitext(file)[1]
	#read the contents of the file into a variable with the same name as the file
	eval(Meta.parse("const $name = readSnaqNetwork(\"$file\")")) #read the whole output
	#put network to variable 'net'
	eval(Meta.parse("net = $name"))
	imagefilename = "$name$r$suffix"
	println(name)
	#Reroot the network using outgroup
	rootatnode!(net, outgroup)
	R"svg"(imagefilename, width=12, height=6)
	R"par"(mar=[0.1,0.1,0.1,0.1])
	plot(net, showgamma=true, tipcex=0.6, tipoffset = 0.2, style = :fulltree);
	R"dev.off()"
end
