using ITensors
using JLD2, FileIO
include("./util.jl")
include("./hamiltonian.jl")
using .util
using .hamiltonian



sweeps = Sweeps(20)
maxdim!(sweeps, 10,20,100,100,200, 500, 1000, 1000, 2000, 2000, 2000, 4000)
cutoff!(sweeps, 1E-30)

N = 16
H_XXZ, sites = Heisenberg_1(N)

e_0, psi_0 = ground_state(H_XXZ, sites, sweeps)
wfs = [psi_0]
e_1, psi_1 = excited_state(wfs, H_XXZ, sites, sweeps)
wfs = [psi_0 psi_1]

s_rel = calculate_relative_entropies(psi_1, psi_0, N)
# println(s_rel)

@save "TFIM_$N.jld2" wfs s_rel
