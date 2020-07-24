using ITensors
using Plots
include("./util.jl")
include("./hamiltonian.jl")
using .util
using .hamiltonian

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
local_dim = dim_dict("$model")

wf1 = util.load_MPS("data/$(model)_$N", "wf1")

plot(calculate_entropies(wf1, N, local_dim), title="$(model)_$N")

savefig("plots/$model $N ground state entropy")
