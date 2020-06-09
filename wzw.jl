using ITensors
# using JLD2, FileIO
using HDF5
include("./util.jl")
include("./hamiltonian.jl")
using .util
using .hamiltonian



sweeps = Sweeps(15)
maxdim!(sweeps, 10,20,100,100,200, 500, 1000, 1000, 2000, 2000, 2000, 4000)
cutoff!(sweeps, 1E-30)

println(ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ()))
N = 16
H_WZW, sites = WZW_2_2(N)

@time e_0, psi_0 = ground_state(H_WZW, sites, sweeps)
wfs = [psi_0]
e_1, psi_1 = excited_state(wfs, H_WZW, sites, sweeps)
wfs = [psi_0 psi_1]

s_rel = calculate_relative_entropies(psi_1, psi_0, N, 3)
println(s_rel)

fname = "data/WZW_$N"
save_MPS(psi_0, fname, "wf1")
save_MPS(psi_1, fname, "wf2")
# save("WZW_$N.jld", "wfs", wfs, "s_rel", s_rel)
