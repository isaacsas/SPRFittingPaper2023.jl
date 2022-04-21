module SPRFitting

using DocStringExtensions
using DataStructures: MutableBinaryHeap, top_with_handle, DataStructures
using DelimitedFiles, Random, StaticArrays
using OnlineStats: Mean, Variance, mean, var, std, nobs, fit!
using UnPack
using JLD, Interpolations, LossFunctions
using BlackBoxOptim: bboptimize, best_candidate, best_fitness
using XLSX, Plots

include("utils.jl")

include("parameters.jl")
export BioPhysParams, SimParams

include("forward_simulator.jl")
export run_spr_sim!

include("callbacks.jl")
export TotalBoundOutputter, TotalAOutputter

include("spr_data.jl")
export AlignedData, get_aligned_data

include("surrogate.jl")
export SurrogateParams, Surrogate, build_surrogate_serial

include("fitting.jl")
export fit_spr_data, bboptpars_to_physpars, visualisefit, savefit

end
