#!/usr/bin/julia

using IceAnalysis

#include("/Users/robinson/models/IceAnalysis.jl/src/IceAnalysis.jl")
#using .IceAnalysis



if length(ARGS) != 1
    throw("Error: path to ensemble must be provided.")
end

path = ARGS[end];


ens = ensemble_check(path);
