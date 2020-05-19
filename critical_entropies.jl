using ITensors
using Plots
using JLD2, FileIO
include("./util.jl")
include("./hamiltonian.jl")
using .util
using .hamiltonian

@load "TFIM_20.jld2" wfs

plot(calculate_entropies(wfs[1], 20), show=true)
