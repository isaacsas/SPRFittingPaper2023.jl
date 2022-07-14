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

function Surrogate(lutdata::AbstractArray; surpars, simpars, surrogate_size=nothing)
    sursize = (surrogate_size === nothing) ? size(lutdata) : surrogate_size
    itp = interpolate(lutdata, BSpline(Linear()))
    @show typeof(itp)
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

######################################## SAVING ##########################################
function save_surrogate_metadata(file, sursize, surpars::SurrogateParams, simpars::SimParams)
    write(file, "surrogate_size", sursize)
    for f in fieldnames(SurrogateParams)
        write(file, string(f), getfield(surpars, f))
    end
    fs = (:N,:tstop,:tstop_AtoB,:tsave,:L)
    for f in fs
        write(file, string(f), getfield(simpars, f))
    end
    write(file, "DIM", getdim(simpars))
    nothing
end

"""
    save_surrogate_metadata(filename, surpars::SurrogateParams, simpars::SimParams)

Create a JLD file with the given surrogate metadata.
"""
function save_surrogate_metadata(filename::String, sursize, surpars::SurrogateParams, simpars::SimParams)
    jldopen(filename, "w") do file
        save_surrogate_metadata(file, sursize, surpars, simpars)
    end
    nothing
end

"""
    save_surrogate(filename, sur::Surrogate)

Create a JLD file with the given surrogate.
"""
function save_surrogate(filename::String, surrogate_size, surpars::SurrogateParams, simpars::SimParams)
    surrogatetable = build_surrogate_serial(surrogate_size, surpars, simpars)
    jldopen(filename, "w") do file
        write(file, "FirstMoment", surrogatetable)
        save_surrogate_metadata(file, surrogate_size, surpars, simpars)
    end
    nothing
end

"""
    save_surrogate_slice(filename, surslice::Matrix, idxstart, idxend)

Create a JLD file with the trajectories (stored in a matrix) corresponding to the parameter
sets `idxstart` through `idxend`.

Notes:
- `surslice` should be number of time samples by `length(idxstart:idxend)`, with each column
  a given solution trajectory for the associated parameter set.
"""
function save_surrogate_slice(filename::String, surmetadata, idxstart, idxend)
    surslice = build_surrogate_slice(surmetadata, idxstart, idxend)

    jldopen(filename, "w") do file
        write(file, "surrogate_slice", surslice)
        write(file, "idxstart", idxstart)
        write(file, "idxend", idxend)
    end
    nothing
end


######################################## BUILDING ##########################################

"""
    build_surrogate_serial(surrogate_size::Tuple, surpars::SurrogateParams, simpars::SimParams;
                           terminator=VarianceTerminator())

Construct the data table array for a new surrogate by varying kon, koff, konb and reach.

Arguments:
- `surrogate_size = ` a Tuple with `(nkon,nkoff,nkonb,nreach)` points to use.
- `surpars = ` the physical parameters to use in the surrogate. It is strongly recommended
  to not change the default values of `antigenconcen`, `antibodyconcen`, or `CP` unless you
  really know what you are doing --  these default values are implicitly assumed in other
  places.
- `simpars = ` the simulation parameters to use (number of particles, simulations, domain
  size, etc).

Keyword Arguments:
- `terminator`, can be used to alter how the number of samples for each parameter set is
  determined. See [`VarianceTerminator`](@ref) for the default values.
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
    surmeans
end

function unpack_metadata(surmetadata)
    surrogate_size = surmetadata["surrogate_size"]
    @unpack logkon_range, logkoff_range, logkonb_range, reach_range = surmetadata
    logkons = collect(range(logkon_range[1], logkon_range[2], length=surrogate_size[1]))
    logkoffs = collect(range(logkoff_range[1], logkoff_range[2], length=surrogate_size[2]))
    logkonbs = collect(range(logkonb_range[1], logkonb_range[2], length=surrogate_size[3]))
    reaches = collect(range(reach_range[1], reach_range[2], length=surrogate_size[4]))
    biopars  = BioPhysParams(; kon=0.0, koff=0.0, konb=0.0, reach=0.0)
    @unpack tstop, tstop_AtoB, tsave = surmetadata
    simpars = SimParams(; tstop, tstop_AtoB, tsave)
    return (; surrogate_size, logkons, logkoffs, logkonbs, reaches, biopars, simpars)
end

function setpars!(biopars, idxs, p)
    konidx,koffidx,konbidx,reachidx = idxs
    biopars.kon  = 10.0 ^ p.logkons[konidx]
    biopars.koff = 10.0 ^ p.logkoffs[koffidx]
    biopars.konb = 10.0 ^ p.logkonbs[konbidx]
    biopars.reach = p.reaches[reachidx]
    nothing
end

"""
    build_surrogate_slice(surmetadata, idxstart, idxend)

Simulates the portion of the surrogate's data corresponding to the linear indices from
idxstart to idxend. Here the indices correspond to the portion of the ordered parameter
space to simulate.
"""
function build_surrogate_slice(surmetadata::Dict, idxstart, idxend; terminator=VarianceTerminator())

    # get the parameters out and processed from the metadata Dict
    p = unpack_metadata(surmetadata)
    @unpack biopars, simpars = p
    @unpack tstop, tstop_AtoB, tsave = simpars

    # output from simulations
    tbo = TotalBoundOutputter(length(tsave))
    idxs = idxstart:idxend
    surrogate_data = zeros(p.surrogate_size[end], length(idxs))
    cidxs = CartesianIndices(p.surrogate_size[1:end-1])

    for (n,idx) in enumerate(idxs)
        # get the parameters for this linear index
        setpars!(biopars, Tuple(cidxs[idx]), p)

        # run a sim and save the average bound
        run_spr_sim!(tbo, biopars, simpars, terminator)
        surrogate_data[:,n] = means(tbo)

        # reset the outputter and terminator
        tbo()
        reset!(terminator)
    end

    surrogate_data
end

function merge_surrogate_slices(p, slicefilebasename, nfiles)
    # get the parameters from the metadata Dict
    @unpack surrogate_size = p

    surmeans = zeros(surrogate_size)
    cidxs = CartesianIndices(surrogate_size[1:end-1])

    for fidx in 1:nfiles
        slicefname = slicefilebasename * "_$fidx.jld"
        data = load(slicefname)
        @unpack surrogate_slice, idxstart, idxend = data

        for (n,idx) in enumerate(idxstart:idxend)
            surmeans[Tuple(cidxs[idx])...,:] = surrogate_slice[:,n]
        end
    end
    surmeans
end

function merge_surrogate_slices(surmetafname::String, slicefilebasename, nfiles; force=false)
    surmetadata = load(surmetafname)
    p = unpack_metadata(surmetadata)
    surmeans = merge_surrogate_slices(p, slicefilebasename, nfiles)

    mergedfname = slicefilebasename * "_merged.jld"
    cp(surmetafname, mergedfname; force)
    jldopen(mergedfname, "r+") do file
        write(file, "FirstMoment", surmeans)
    end
    return mergedfname
end