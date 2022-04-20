"""
$(TYPEDEF)

Ranges of parameters used in the surrogate.

# Fields
$(FIELDS)
"""
Base.@kwdef struct SurrogateRanges
    """log₁₀ space range of `kon`."""
    logkon_range::Tuple{Float64,Float64}
    """log₁₀ space range of `koff`."""
    logkoff_range::Tuple{Float64,Float64}
    """log₁₀ space range of `konb`."""
    logkonb_range::Tuple{Float64,Float64}
    """linear range of reach values"""
    reach_range::Tuple{Float64,Float64}
end

"""
$(TYPEDEF)

Biophysical parameters used in the surrogate.

# Fields
$(FIELDS)
"""
Base.@kwdef struct SurrogateParams
    """log₁₀ space `kon`s."""
    logkonvec::Vector{Float64}
    """log₁₀ space `koff`s."""
    logkoffvec::Vector{Float64}
    """log₁₀ space `konb`s."""
    logkonbvec::Vector{Float64}
    """linear space `reach`s"""
    reachvec::Vector{Float64}
    """Surrogate's internal antigen concentration in μM (default is `DEFAULT_SIM_ANTIGENCONCEN`)."""
    antigenconcen = DEFAULT_SIM_ANTIGENCONCEN
    """Surrogate's inteneral antibody concentration in μM (default is 1.0)."""
    antibodyconcen = 1.0
    """Surrogate's internal CP value (default is 1.0)"""
    CP = 1.0
end

"""
SurrogateParams(param_ranges::SurrogateRanges; numkon, numkoff, numkonb, numreach, kwargs...)
"""
function SurrogateParams(param_ranges::SurrogateRanges; numkon, numkoff, numkonb, numreach, kwargs...)
    r = param_ranges.logkon_range
    logkonvec = range(r[1], r[2], length=numkon)
    r = param_ranges.logkoff_range
    logkoffvec = range(r[1], r[2], length=numkoff)
    r = param_ranges.logkonb_range
    logkonbvec = range(r[1], r[2], length=numkonb)
    r = param_ranges.reach_range
    reachvec = range(r[1], r[2], length=numreach)

    SurrogateParams(; logkonvec, logkoffvec, logkonbvec, reachvec, kwargs...)
end
"""
$(TYPEDEF)

Surrogate parameters and data.

# Fields
$(FIELDS)

Keyword Arguments:
- `param_ranges = ` a [`SurrogateRanges`](@ref) defining the parameter ranges.
- `lutfile = ` the name of the file storing a surrogate to load.
- `lutdata = ` `AbstractArray` representing the raw data points to build the
  surrogate from.
- `antigenconcen = DEFAULT_SIM_ANTIGENCONCEN` the surrogate's internal antigen
  concentration.
"""
Base.@kwdef mutable struct Surrogate{S,T,U,V}
    """[`SurrogateRanges`](@ref) representing ranges for each parameter."""
    param_ranges::SurrogateRanges
    """Number of points for each coordinate within the surrogate lookup table."""
    surrogate_size::S
    """Interpolation of the surrogate lookup table."""
    itp::T
    """Biophysical parameters used in the surrogate."""
    biopars::SurrogateParams
    """Simulation parameters used in the surrogate."""
    simpars::SimParams{U,V}
end

function Surrogate(lutdata::AbstractArray; kwargs...)
    surrogate_size = size(lutdata)
    itp = interpolate(lutdata, BSpline(Linear()))
    Surrogate(; surrogate_size, itp, kwargs...)
end

function Surrogate(lutfile::String; rungc=true, kwargs...)
    lutdata = load(lutfile)["FirstMoment"]
    surrogate = Surrogate(lutdata; kwargs...)

    # release the lookup table from memory
    if rungc
        lutdata = nothing
        GC.gc()
    end

    surrogate
end


# function build_surrogate(pars::SurrogateParams, biopars, simpars)

# end

