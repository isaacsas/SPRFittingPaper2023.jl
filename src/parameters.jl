# This is the default internal antigen concentration 
# used in the surrogates and simulations.
# 125.23622683286348 μM
const DEFAULT_SIM_ANTIGENCONCEN = 500/149/26795*1000000  


"""
$(TYPEDEF)

Biophysical parameters to use in forward simulations.

# Fields
$(FIELDS)
"""
Base.@kwdef mutable struct BioPhysParams{T <: Number}
    """A --> B rate, units of (μM s)⁻¹"""
    kon::T
    """B --> A rate, units of s⁻¹"""
    koff::T
    """A+B --> C rate, units of s⁻¹"""
    konb::T
    """Reach of A+B --> C reaction, units of nm"""
    reach::T
    """CP factor (default = 1.0)"""
    CP::T = 1.0
    """Concentration of antigen (default = `DEFAULT_SIM_ANTIGENCONCEN` μM)"""
    antigenconcen::T = DEFAULT_SIM_ANTIGENCONCEN
    """Concentration of antibodies, should be consistent with kon's units (default = 1.0 μM)"""
    antibodyconcen::T = 1.0
end

"""
    biopars_from_fitting_vec(p; antibodyconcen=1.0, antigenconcen=1.0)

Generate a [`BioPhysParams`](@ref) struct from a parameter vector from fitting.

Notes: 
- Assumed that 
    ``
    p = [\\log_{10}(\\text{kon}), \\log_{10}(\\text{koff}), \\log_{10}(\\text{kon}_{\\text{b}}), \\text{reach}, \\log_{10}(\\text{CP})] 
    ``
"""
function biopars_from_fitting_vec(p; antibodyconcen=1.0, antigenconcen=DEFAULT_SIM_ANTIGENCONCEN)
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


"""
$(TYPEDEF)

Biophysical parameters to use in forward simulations.

# Fields
$(FIELDS)

Keyword Arguments:
- All fields have a corresponding kwarg.
- `antigenconcen = DEFAULT_SIM_ANTIGENCONCEN`, the default antigenconcen in
  units of (nm)⁻³ (unless using `convert_agc_units=false`).
- `DIM = 3` the dimension of the underlying space (i.e. 2 or 3). 
- `convert_agc_units = true`, set to false to disable the conversion of the
  antigen units from assumed units of μM to (nm)⁻³, see below. In this case
  ``L`` is calculated using the antigen concentration value directly. 
- `resample_initlocs = true`, if set to false `initlocs` will be constant,
  simply reusing the initial value sampled in the `SimParams` constructor.

Notes:
- Uses the antigen concentration to determine ``L``, and so needs to be updated
  if this concentration changes. The antigen concentration is first converted
  from assumed units of  μM to units of (nm)⁻³. Then 
  ```math
    L = \\left(\\frac{N}{[\\text{antigen}] \\text{ in (nm)}^{-3}}\\right)^{\\tfrac{1}{DIM}}
  ```
"""
mutable struct SimParams{T<:Number,DIM}
    """Number of particles (default `= 1000`)"""
    N::Int
    """Time to end simulations at (default `= 600.0`)"""
    tstop::T
    """Time to turn off A --> B reaction, i.e. time the antibody bath is removed. (default `= Inf`)"""
    tstop_AtoB::T
    """Times to save data at (default = `nothing`)"""
    tsave::Vector{T}
    """Domain Length"""
    L::T 
    """Initial Particle Positions (default = uniform in ``[-L/2,L/2]^{\\text{DIM}}``)"""
    initlocs::Vector{SVector{DIM,Float64}}
    """Resample `initlocs` every simulation (default `= true`)"""
    resample_initlocs::Bool
    """Number of simulations in runsim (to control sampling error) (default `= 1000`)"""
    nsims::Int
end

function Base.show(io::IO, ::MIME"text/plain", sps::SimParams)   
    @unpack N,tstop,tstop_AtoB,tsave,L,resample_initlocs,nsims = sps 
    println(io, summary(sps))
    println(io, "number of particles (N) = ", N)
    println(io, "tstop = ", tstop)
    println(io, "tstop_AtoB = ", tstop_AtoB)
    println(io, "number of save points = ", length(tsave))
    println(io, "domain length (L) = ", L)
    println(io, "resample_initlocs = ", resample_initlocs)
    print(io, "nsims = ", nsims)
end

# assumes antigenconcen is μM
function SimParams(; antigenconcen=DEFAULT_SIM_ANTIGENCONCEN,
                    N=1000, 
                    tstop=600.0, 
                    tstop_AtoB=Inf, 
                    dt=1.0, 
                    tsave=nothing,
                    L=nothing, 
                    DIM=3,
                    initlocs=nothing, 
                    resample_initlocs=true,
                    nsims=1000,
                    convert_agc_units=true)

    # convert agc to (nm)⁻³                
    agc = convert_agc_units ? muM_to_inv_cubic_nm(antigenconcen) : antigenconcen
    Lv  = isnothing(L) ? (N / agc)^(1/DIM) : L
    initlocsv = if (initlocs === nothing) 
        [(Lv .* rand(SVector{DIM,Float64}) .- Lv/2) for _ in 1:N] 
    else 
         initlocs    
    end

    # if saving times not given use dt to determine
    tsavev = (tsave === nothing) ? collect(range(0.0, tstop, step=dt)) : tsave

    SimParams{typeof(tstop),DIM}(N,tstop,tstop_AtoB,tsavev,Lv,initlocsv,resample_initlocs,nsims)
end

function SimParams(biopars::BioPhysParams; kwargs...)
    SimParams(biopars.antigenconcen; kwargs...)
end

######################## helpful accessors #####################

getdim(simpars::SimParams{T,DIM}) where {T<:Number,DIM} = DIM

getantigenconcen(simpars) = simpars.N / (simpars.L^getdim(simpars))
