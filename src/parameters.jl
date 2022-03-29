
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


# Note the following depends on antigen concentration to set L, so needs to be recalculated if it changes!
mutable struct SimParams{T<:Number,U,V<:AbstractArray{T}}
    """Number of particles"""
    N::Int
    """Time to end simulations at"""
    tstop::T
    """Time to turn off A --> B reaction"""
    ts_1::T
    """How often to save data in time"""
    dt::T
    """Time vector for interpolation"""
    t::U 
    """Domain Length"""
    L::T 
    """Initial Particle Positions"""
    initlocs::V
    """Resample initlocs every simulation"""
    resample_initlocs::Bool
    """Number of simulations in runsim (to control sampling error)"""
    nsims::Int
end

# converts antigen concentration to number of units
function SimParams(antigenconcen; 
                   N=1000, 
                   tstop=600.0, 
                   ts_1=150.0, 
                   dt=1.0, 
                   t=nothing, 
                   L=nothing, 
                   initlocs=nothing, 
                   resample_initlocs=true,
                   nsims=1000)
    tv  = isnothing(t) ? range(1.0, tstop+1, step=1.0) : t    
    Lv  = isnothing(L) ? sqrt(N / (antigenconcen)) : L
    initlocsv = isnothing(initlocs) ? (Lv .* rand(2,N) .- Lv/2) : initlocs    
    SimParams(N,tstop,ts_1,dt,tv,Lv,initlocsv,resample_initlocs,nsims)
end

function SimParams(biopars::BioPhysParams; kwargs...)
    SimParams(biopars.antigenconcen; kwargs...)
end

#################################