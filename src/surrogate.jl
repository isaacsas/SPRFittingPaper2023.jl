"""
$(TYPEDEF)

Biophysical parameters used in the surrogate.

# Fields
$(FIELDS)
"""
Base.@kwdef struct SurrogateParams{T <: Number}
    """log₁₀ space range of `kon`, first order rate in surrogate (units of s⁻¹)."""
    logkon_range::Tuple{T,T}
    """log₁₀ space range of `koff` (units of s⁻¹)."""
    logkoff_range::Tuple{T,T}
    """log₁₀ space range of `konb`, Doi association rate in surrogate (units of s⁻¹)."""
    logkonb_range::Tuple{T,T}
    """linear range of reach values (units of nm)"""
    reach_range::Tuple{T,T}
    """Surrogate's internal antigen concentration in μM (default is `DEFAULT_SIM_ANTIGENCONCEN`)."""
    antigenconcen::T = DEFAULT_SIM_ANTIGENCONCEN
    """Surrogate's internal antibody concentration in μM (default is 1.0)."""
    antibodyconcen::T = 1.0
    """Surrogate's internal CP value (default is 1.0)"""
    CP::T = 1.0
end

function Base.show(io::IO,  ::MIME"text/plain", sur::SurrogateParams)   
    println(io, summary(sur))
    println(io, "logkon_range = ", sur.logkon_range)
    println(io, "logkoff_range = ", sur.logkoff_range)
    println(io, "logkonb_range = ", sur.logkonb_range)
    println(io, "reach_range = ", sur.reach_range)
    println(io, "[antigen] = ", sur.antigenconcen)
    println(io, "[antibody] = ", sur.antibodyconcen)
    print(io, "CP = ", sur.CP)    
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

"""
    save_surrogate(filename, sur::Surrogate)

Create a JLD file with the given surrogate.
"""
function save_surrogate(filename, sur::Surrogate{S,T}) where {S, T <: Array}

    jldopen(filename, "w") do file
        write(file, "FirstMoment", sur.itp)
    
        for f in fieldnames(SurrogateParams)
            write(file, string(f), getfield(sur.surpars, f))
        end

        fs = (:N,:tstop,:tstop_AtoB,:tsave,:L)
        for f in fs
            write(file, string(f), getfield(sur.simpars, f))
        end
        write(file, "DIM", getdim(sur.simpars))
    end
    
    nothing
end

"""
    build_surrogate_serial(surrogate_size::Tuple, surpars::SurrogateParams, simpars::SimParams; 
                           terminator=VarianceTerminator())

Creates a new surrogate varying kon, koff, konb and reach uniformly in log
space.

Arguments:
- `surrogate_size = ` a Tuple with `(nkon,nkoff,nkonb,nreach)` points to use.
- `surpars = ` the physical parameters to use in the surrogate. It is strongly
  recommended to not change the default values of `antigenconcen`,
  `antibodyconcen`, or `CP` unless you really know what you are doing --  these
  default values are implicitly assumed in other places.
- `simpars = ` the simulation parameters to use (number of particles,
  simulations, domain size, etc).

Keyword Arguments:
- `terminator`, can be used to alter how the number of samples for each
  parameter set is determined. See [`VarianceTerminator`](@ref) for the default
  values.
"""
function build_surrogate_serial(surrogate_size::Tuple, surpars::SurrogateParams, simpars::SimParams; 
                                terminator=VarianceTerminator())
    @assert surrogate_size[end] == length(simpars.tsave)

    # reach of parameters we vary
    logkons  = range(surpars.logkon_range[1], surpars.logkon_range[2], length=surrogate_size[1])
    logkoffs = range(surpars.logkoff_range[1], surpars.logkoff_range[2], length=surrogate_size[2])
    logkonbs = range(surpars.logkonb_range[1], surpars.logkonb_range[2], length=surrogate_size[3])
    reachs   = range(surpars.reach_range[1], surpars.reach_range[2], length=surrogate_size[4])    
    biopars  = BioPhysParams(; kon=0.0, koff=0.0, konb=0.0, reach=0.0)
    
    # output from simulations
    tbo      = TotalBoundOutputter(length(simpars.tsave))
    surmeans = zeros(surrogate_size)
    
    # run and save the simulation results 
    for (i4,reach) in enumerate(reachs)
        biopars.reach = reach 

        for (i3,logkonb) in enumerate(logkonbs)
            biopars.konb = 10.0^logkonb

            for (i2,logkoff) in enumerate(logkoffs)
                biopars.koff = 10.0^logkoff

                for (i1,logkon) in enumerate(logkons)
                    biopars.kon = 10.0^logkon
                    
                    run_spr_sim!(tbo, biopars, simpars, terminator)
                    means!(view(surmeans,i1,i2,i3,i4,:), tbo)

                    # reset outputter and terminator
                    tbo()                                    
                    reset!(terminator)
                end
            end
        end
    end

    Surrogate(surrogate_size, surmeans, surpars, simpars)
end
