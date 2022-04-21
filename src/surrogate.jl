"""
$(TYPEDEF)

Biophysical parameters used in the surrogate.

# Fields
$(FIELDS)
"""
Base.@kwdef struct SurrogateParams{T <: Number}
    """log₁₀ space range of `kon`."""
    logkon_range::Tuple{T,T}
    """log₁₀ space range of `koff`."""
    logkoff_range::Tuple{T,T}
    """log₁₀ space range of `konb`."""
    logkonb_range::Tuple{T,T}
    """linear range of reach values"""
    reach_range::Tuple{T,T}
    """Surrogate's internal antigen concentration in μM (default is `DEFAULT_SIM_ANTIGENCONCEN`)."""
    antigenconcen::T = DEFAULT_SIM_ANTIGENCONCEN
    """Surrogate's intenral antibody concentration in μM (default is 1.0)."""
    antibodyconcen::T = 1.0
    """Surrogate's internal CP value (default is 1.0)"""
    CP::T = 1.0
end


"""
$(TYPEDEF)

Surrogate parameters and data.

# Fields
$(FIELDS)

Arguments (one of the following two):
- `lutfile = ` the name of the file storing a surrogate to load.
- `lutdata = ` `AbstractArray` representing the raw data points to build the
  surrogate from.

Notes:
- All fields can also be passed as a keyword arg.
"""
Base.@kwdef mutable struct Surrogate{S,T,U,V,W}
    """Number of points for each coordinate within the surrogate lookup table."""
    surrogate_size::S
    """Interpolation of the surrogate lookup table."""
    itp::T
    """[`SurrogateParams`](@ref) used in the surrogate."""
    surpars::SurrogateParams{U}
    """Simulation parameters used in the surrogate."""
    simpars::SimParams{V,W}
end

function Surrogate(lutdata::AbstractArray; surpars, simpars)
    surrogate_size = size(lutdata)
    itp = interpolate(lutdata, BSpline(Linear()))
    Surrogate(; surrogate_size, itp, surpars, simpars)
end

function Surrogate(lutfile::String; rungc=true, surpars=nothing, simpars=nothing)
    lutdata = load(lutfile)
    itpdata = lutdata["FirstMoment"]

    if surpars === nothing
        fs = fieldnames(SurrogateParams)
        all(f -> haskey(lutdata, string(f)), fs) || error("Not all SurrogateParams fields are in the surrogate file.")
        kwargs   = NamedTuple(f => lutdata[string(f)] for f in fs)
        surpars′ = SurrogateParams(; kwargs...)
    else
        surpars′ = surpars
    end

    if simpars === nothing
        fs = (:N,:tstop,:tstop_AtoB,:tsave,:L,:DIM)
        all(f -> haskey(lutdata, string(f)), fs) || error("Not all needed SimParams fields are in the surrogate file.")
        kwargs   = NamedTuple(f => lutdata[string(f)] for f in fs)
        simpars′ = SimParams(; antigenconcen=surpars′.antigenconcen, kwargs...)
    else
        simpars′ = simpars
    end

    surrogate = Surrogate(itpdata; surpars=surpars′, simpars=simpars′)

    # release the lookup table from memory
    if rungc
        lutdata = nothing
        GC.gc()
    end

    surrogate
end


# function build_surrogate(pars::SurrogateParams, biopars, simpars)

# end

