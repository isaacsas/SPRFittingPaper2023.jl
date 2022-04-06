
Base.@kwdef mutable struct BioPhysParams{T <: Number}
    """A --> B rate, units of concentration per time"""
    kon::T
    """B --> A rate, units of per time"""
    koff::T
    """A+B --> C rate, units of per time"""
    konb::T
    """Reach of A+B --> C reaction"""
    reach::T
    """CP factor"""
    CP::T = 1.0
    """Concentration of antigen"""
    antigenconcen::T = 1.0
    """Concentration of antibodies"""
    antibodyconcen::T = 1.0
end

# generate BioPhysParams from a parameter vector from fitting
# assumes: 
# p = [log10(kon), log10(koff), log10(konb), reach, log10(CP)] 
function biopars_from_fitting_vec(p; antibodyconcen=1.0, antigenconcen=1.0)
    kon   = 10^p[1]
    koff  = 10^p[2]
    konb  = 10^p[3]
    reach = p[4]
    CP    = 10^p[5]
    BioPhysParams(kon,koff,konb,reach,CP,antigenconcen,antibodyconcen)
end

function Base.show(io::IO,  ::MIME"text/plain", bps::BioPhysParams)   
    @unpack kon,koff,konb,reach,CP,antigenconcen,antibodyconcen = bps 
    println(io, summary(bps))
    println(io, "kon = ", kon)
    println(io, "koff = ", koff)
    println(io, "konb = ", konb)
    println(io, "reach = ", reach)
    println(io, "CP = ", CP)
    println(io, "[antigen] = ", antigenconcen)
    print(io, "[antibody] = ", antibodyconcen)
end

# Note the following depends on antigen concentration to set L, so needs to be recalculated if it changes!
mutable struct SimParams{T<:Number,DIM}
    """Number of particles"""
    N::Int
    """Time to end simulations at"""
    tstop::T
    """Time to turn off A --> B reaction"""
    tstop_AtoB::T
    """How often to save data in time"""
    dt::T
    """Domain Length"""
    L::T 
    """Initial Particle Positions"""
    initlocs::Vector{SVector{DIM,Float64}}
    """Resample initlocs every simulation"""
    resample_initlocs::Bool
    """Number of simulations in runsim (to control sampling error)"""
    nsims::Int
end

function Base.show(io::IO, ::MIME"text/plain", sps::SimParams)   
    @unpack N,tstop,tstop_AtoB,dt,L,resample_initlocs,nsims = sps 
    println(io, summary(sps))
    println(io, "number of particles (N) = ", N)
    println(io, "tstop = ", tstop)
    println(io, "tstop_AtoB = ", tstop_AtoB)
    println(io, "save frequency (dt) = ", dt)
    println(io, "domain length (L) = ", L)
    println(io, "resample_initlocs = ", resample_initlocs)
    print(io, "nsims = ", nsims)
end

# converts antigen concentration to number of units
function SimParams(antigenconcen; 
                   N=1000, 
                   tstop=600.0, 
                   tstop_AtoB=Inf, 
                   dt=1.0, 
                   L=nothing, 
                   DIM=3,
                   initlocs=nothing, 
                   resample_initlocs=true,
                   nsims=1000)
    Lv  = isnothing(L) ? (N / antigenconcen)^(1/DIM) : L
    initlocsv = if (initlocs === nothing) 
        [(Lv .* rand(SVector{DIM,Float64}) .- Lv/2) for _ in 1:N] 
    else 
         initlocs    
    end
    SimParams{typeof(tstop),DIM}(N,tstop,tstop_AtoB,dt,Lv,initlocsv,resample_initlocs,nsims)
end

function SimParams(biopars::BioPhysParams; kwargs...)
    SimParams(biopars.antigenconcen; kwargs...)
end

