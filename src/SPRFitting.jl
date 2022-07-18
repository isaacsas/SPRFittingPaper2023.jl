module SPRFitting

using DocStringExtensions, UnPack
using DataStructures: MutableBinaryHeap, top_with_handle, DataStructures
using CSV: File, CSV
using Tables, Random, StaticArrays
using OnlineStats: Mean, Variance, mean, var, std, nobs, fit!
using JLD, Interpolations
using BlackBoxOptim: bboptimize, best_candidate, best_fitness
using XLSX, Plots
using Optimization

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
export SurrogateParams, Surrogate, save_surrogate, save_surrogate_metadata,
       save_surrogate_slice, merge_surrogate_slices

include("fitting.jl")
export fit_spr_data, bboptpars_to_physpars, visualisefit, savefit

include("monovalent_fitting.jl")
export monovalent_fit_spr_data

end
