################ STRUCTURES FOR SAVING OUTPUT ################

"""
$(TYPEDEF)

Callback that saves the amount of bound antibodies (i.e. monovalently +
bivalently bound) at each save time.

# Fields
$(FIELDS)
"""
struct TotalBoundOutputter{T}
    bindcnt::Vector{T}
end

# construct for `N` time points 
function TotalBoundOutputter(N) 
    v = [Variance() for _ in 1:N]
    TotalBoundOutputter(v)
end

# this just sums up the amount of B and C at each time, and saves the mean/variance
@inline function (o::TotalBoundOutputter)(tsave, copynumbers, biopars, simpars)
    @unpack CP = biopars
    @unpack N  = simpars

    @inbounds for i in axes(copynumbers,2)
        totalbound = (CP/N) * (copynumbers[2,i] + copynumbers[3,i])
        fit!(o.bindcnt[i], totalbound)
    end

    nothing
end

# called after all simulations finish
@inline function (o::TotalBoundOutputter)(biopars, simpars)
    nothing 
end

# this is used to reset the output object prior to a new simulation (i.e. with
# new parameters)
@inline function (o::TotalBoundOutputter)()
    for i in eachindex(o.bindcnt)
        o.bindcnt[i] = Variance()
    end
    nothing 
end

"""
    means(o::TotalBoundOutputter)
    
Vector of means.
"""
means(o::TotalBoundOutputter) = mean.(o.bindcnt)

"""
    means!(m, o::TotalBoundOutputter)
    
In-place vector of means.
"""
function means!(m, o::TotalBoundOutputter) 
    for i in eachindex(m)
        m[i] = mean(o.bindcnt[i])
    end
    nothing
end

"""
    vars(o::TotalBoundOutputter)
    
Vector of variances.
"""
vars(o::TotalBoundOutputter) = var.(o.bindcnt)

"""
    means(o::TotalBoundOutputter)
    
Vector of standard errors.
"""
sems(o::TotalBoundOutputter) = [std(m)/sqrt(nobs(m)) for m in o.bindcnt]


"""
$(TYPEDEF)

Callback that saves the amount of unbound antigen at each save time.

# Fields
$(FIELDS)
"""
struct TotalAOutputter{T}
    bindcnt::Vector{T}
end

# create the Outputter via knowing N, where 
# N = number of time points to save at
function TotalAOutputter(N::Int) 
    m = [Mean() for _ in 1:N]
    TotalAOutputter(m)
end

# this is called after each individual stochastic simulation to handle
# processing the counts of each species at each time. here we just cumulatively
# sum up the amount of A at each time
@inline function (o::TotalAOutputter)(tsave, copynumbers, biopars, simpars)    
    @unpack CP = biopars
    @unpack N = simpars
    @unpack bindcnt = o

    for i in eachindex(bindcnt)
        fit!(bindcnt[i], (CP/N) * copynumbers[1,i])
    end 
    nothing
end

# this is called once all simulations finish, for a given parameter set, to
# finalize the output. Here we rescale to get a normalized average fraction of
# particles in the A state.
@inline function (o::TotalAOutputter)(biopars, simpars)
    nothing
end

# this is used to reset the output object prior to a new simulations (i.e. with
# new parameters)
@inline function (o::TotalAOutputter)()
    for i in eachindex(o.bindcnt)
        o.bindcnt[i] = Mean()
    end
    nothing 
end

"""
    means(o::TotalAOutputter)
    
Vector of means.
"""
means(o::TotalAOutputter) = mean.(o.bindcnt)

#################################################################

# Callbacks for determining if need more simulations or not

"""
$(TYPEDEF)

Callback that stops simulating when the desired number of simulations is
reached.

# Fields
$(FIELDS)
"""
mutable struct SimNumberTerminator
    """How many simulations have been completed"""
    num_completed_sims::Int
end

SimNumberTerminator() = SimNumberTerminator(0)

# called before a simulation to see if it should be executed
@inline function isnotdone(nt::SimNumberTerminator, biopars, simpars)
    nt.num_completed_sims < simpars.nsims
end

# called after a simulation to update the SimNumberTerminator
@inline function update!(f::SimNumberTerminator, outputter, biopars, simpars)
    f.num_completed_sims += 1
    nothing 
end

# called to reset the SimNumberTerminator
@inline function reset!(f::SimNumberTerminator)
    f.num_completed_sims = 0
    nothing 
end


"""
$(TYPEDEF)

Callback that stops simulating when the variance in bound antibodies becomes
sufficiently small or a maximum number of simulations is reached.

# Fields
$(FIELDS)
"""
Base.@kwdef mutable struct VarianceTerminator
    """The tolerance below which to stop simulating (default = .01)."""
    ssetol::Float64 = 0.01
    """True if should keep iterating."""    
    notdone::Bool = true
    """Minimum number of sims to run (default = 15)."""
    minsims::Int = 15
    """Maximum number of sims to run (default = 250)."""
    maxsims::Int = 250
end


# called before a simulation to see if it should be executed
@inline function isnotdone(vt::VarianceTerminator, biopars, simpars)
    vt.notdone
end

# called after a simulation to update the VarianceTerminator
@inline function update!(vt::VarianceTerminator, outputter, biopars, simpars)
    obs   = outputter.bindcnt
    nsims = nobs(obs[1])

    vt.notdone = if nsims < vt.minsims
        true
    elseif nsims >= vt.maxsims
        false
    else
        any(std(o) > vt.ssetol*mean(o)*sqrt(nobs(o)) for o in obs)
    end    

    nothing 
end

# called to reset the VarianceTerminator
@inline function reset!(vt::VarianceTerminator)
    vt.notdone = true
    nothing 
end
