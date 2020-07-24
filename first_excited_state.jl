using ITensors
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
    "--model"
      help = "Name of spin chain"
      arg_type = String
      default = "TFIM"
end

args = parse_args(s)
N = args["N_sites"]
model = args["model"]


sweeps = Sweeps(20)
maxdim!(sweeps, 10,20,100,100,200, 500, 1000, 1000, 2000, 2000, 2000, 4000)
cutoff!(sweeps, 1E-50)

H, sites = H_dict[model](N)

e_0, psi_0 = ground_state(H, sites, sweeps)
wfs = [psi_0]
e_1, psi_1 = excited_state(wfs, H, sites, sweeps)
wfs = [psi_0 psi_1]
println("Success!")
println(wfs)



fname = "$model_$N"
save_MPS(wfs[1], fname, "wf1")
save_MPS(wfs[2], fname, "wf2")
