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
end

function Surrogate(param_ranges::SurrogateRanges, LUTdata::AbstractArray)
    surrogatesize = size(LUTdata)
    itp = interpolate(LUTdata, BSpline(Linear()))
    Surrogate{typeof(surrogatesize),typeof(itp)}(param_ranges, surrogatesize, itp)
end

function Surrogate(param_ranges::SurrogateRanges, lutfile::String; rungc=true)
    lutdata = load(lutfile)["FirstMoment"]    
    surrogate = Surrogate(param_ranges, lutdata)    

    # release the lookup table from memory
    if rungc
        lutdata = nothing
        GC.gc()
    end

    surrogate
end