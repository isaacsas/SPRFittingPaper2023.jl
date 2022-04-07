Base.@kwdef struct SurrogateRanges
    logkon_range::Tuple{Float64,Float64}
    logkoff_range::Tuple{Float64,Float64}
    logkonb_range::Tuple{Float64,Float64}
    """Linear range of reach values"""
    reach_range::Tuple{Float64,Float64}
    logCP_range::Tuple{Float64,Float64}
end

Base.@kwdef mutable struct Surrogate{S,T} 
    """Log space ranges for each parameter (lin space for reach)"""
    param_ranges::SurrogateRanges
    """Number of points for each coordinate within the surrogate lookup table"""
    surrogate_size::S
    """Interpolation of the surrogate lookup table"""
    itp::T
    """Surrogate's internal antigen concentration in μM"""
    antigenconcen
end

function Surrogate(param_ranges::SurrogateRanges, LUTdata::AbstractArray;
                   antigenconcen=DEFAULT_SIM_ANTIGENCONCEN)
    surrogatesize = size(LUTdata)
    itp = interpolate(LUTdata, BSpline(Linear()))
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