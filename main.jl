using ITensors
using JLD2, FileIO
include("./util.jl")
include("./hamiltonian.jl")
using .util
using .hamiltonian
using ArgParse


s = ArgParseSettings()
@add_arg_table! s begin
    "--N_sites", "-N"
    help = "Number of sites"
    arg_type = Int
    default = 8
end

args = parse_args(s)
N = args["N_sites"]


assert false

sweeps = Sweeps(20)
maxdim!(sweeps, 10,20,100,100,200, 500, 1000, 1000, 2000, 2000, 2000, 4000)
cutoff!(sweeps, 1E-50)

H_TFIM, sites = TFIM_Hamiltonian(N)

e_0, psi_0 = ground_state(H_TFIM, sites, sweeps)
wfs = [psi_0]
e_1, psi_1 = excited_state(wfs, H_TFIM, sites, sweeps)
wfs = [psi_0 psi_1]

s_rel = calculate_relative_entropies(psi_1, psi_0, N)
# println(s_rel)

fname = "TFIM_$N"
save_MPS(wfs[1], fname, "wf1")
save_MPS(wfs[2], fname, "wf2")
end
