using ITensors
using Plots
# using JLD2, FileIO
using JLD
include("./util.jl")
include("./hamiltonian.jl")
using .util
using .hamiltonian

wfs = load("WZW_16.jld", "wfs")

plot(calculate_entropies(wfs[1], 20), title="WZW_16")

savefig("../plots/WZW 16 plot")
