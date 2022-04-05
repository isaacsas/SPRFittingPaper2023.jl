module SPRFitting

using DataStructures, Random, StaticArrays
using UnPack


include("parameters.jl")
include("forward_simulator.jl")

export BioPhysParams, SimParams, run_spr_sim!

end
