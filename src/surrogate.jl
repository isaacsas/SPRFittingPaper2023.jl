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
Base.@kwdef mutable struct Surrogate{S,T} 
    """[`SurrogateRanges`](@ref) representing ranges for each parameter."""
    param_ranges::SurrogateRanges
    """Number of points for each coordinate within the surrogate lookup table."""
    surrogate_size::S
    """Interpolation of the surrogate lookup table."""
    itp::T
    """Surrogate's internal antigen concentration in μM."""
    antigenconcen
end

function Surrogate(param_ranges::SurrogateRanges, lutdata::AbstractArray;
                   antigenconcen=DEFAULT_SIM_ANTIGENCONCEN)
    surrogatesize = size(lutdata)
    itp = interpolate(lutdata, BSpline(Linear()))
    Surrogate{typeof(surrogatesize),typeof(itp)}(param_ranges, surrogatesize, itp, antigenconcen)
end

function Surrogate(param_ranges::SurrogateRanges, lutfile::String; 
                   antigenconcen=DEFAULT_SIM_ANTIGENCONCEN, rungc=true)
    lutdata = load(lutfile)["FirstMoment"]    
    surrogate = Surrogate(param_ranges, lutdata; antigenconcen)    

    # release the lookup table from memory
    if rungc
        lutdata = nothing
        GC.gc()
    end

    surrogate
end

