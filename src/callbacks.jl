# structures for saving outputs:

################ SAVING TOTAL BOUND AT EACH TIME ################
"""
$(TYPEDEF)

Callback that saves the amount of bound antibodies (i.e. monovalently +
bivalently bound) at each save time.

# Fields
$(FIELDS)
"""
struct TotalBoundOutputter
    """Amount bound"""
    bindcnt::Vector{Float64}
    """Scaled standard error in amount bound"""
    bindcntsse::Vector{Float64}
end

# allows to create
function TotalBoundOutputter(N::Int) 
    TotalBoundOutputter(zeros(N),zeros(N))
end

# this just sums up the amount of B and C at each time
@inline function (o::TotalBoundOutputter)(tsave, copynumbers, biopars, numpars)
    @unpack CP = biopars
    @unpack N  = numpars

    @inbounds for i in axes(copynumbers,2)
        totalbound       = (CP/N) * (copynumbers[2,i] + copynumbers[3,i])
        o.bindcnt[i]    += totalbound
        o.bindcntsse[i] += totalbound * totalbound
    end
    nothing
end

# this is called once all simulations finish to finalize the output
@inline function (o::TotalBoundOutputter)(biopars, numpars)
    @unpack nsims = numpars

    # mean 
    o.bindcnt ./= nsims

    # variance
    @. o.bindcntsse = abs(o.bindcntsse/(nsims-1) - (nsims/(nsims-1))*(o.bindcnt^2))

    # scaled standard error 
    @. o.bindcntsse = sqrt(o.bindcntsse/nsims) / o.bindcnt
    for i in eachindex(o.bindcntsse)
        (o.bindcnt[i] < eps()) && (o.bindcntsse[i] = 0.0)
    end

    nothing
end

# this is used to reset the output object prior to a new simulations (i.e. with
# new parameters)
@inline function (o::TotalBoundOutputter)()
    o.bindcnt .= 0
    o.bindcntsse .= 0
    nothing 
end


################################################################

############ SAVING TOTAL AMOUNT OF A AT EACH TIME #############

"""
$(TYPEDEF)

Callback that saves the amount of unbound antigen at each save time.

# Fields
$(FIELDS)
"""
struct TotalAOutputter
    bindcnt::Vector{Float64}
end

# create the Outputter via knowing N, where 
# N = number of time points to save at
function TotalAOutputter(N::Int) 
    TotalAOutputter(zeros(N))
end

# this is called after each individual stochastic simulation to handle
# processing the counts of each species at each time. here we just cumulatively
# sum up the amount of A at each time
@inline function (o::TotalAOutputter)(tsave, copynumbers, biopars, numpars)    
    o.bindcnt .+= @view copynumbers[1,:]
    nothing
end

# this is called once all simulations finish, for a given parameter set, to
# finalize the output. Here we rescale to get a normalized average fraction of 
# particles in the A state.
@inline function (o::TotalAOutputter)(biopars, numpars)
    @unpack CP = biopars
    @unpack N,nsims = numpars
    o.bindcnt .*= (CP/(N*nsims))
    nothing
end

# this is used to reset the output object prior to a new simulations (i.e. with
# new parameters)
@inline function (o::TotalAOutputter)()
    o.bindcnt .= 0
    nothing 
end

#################################################################

# callbacks for determining if need more simulations or not

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
@inline function isnotdone(nt::SimNumberTerminator, biopars, numpars)
    nt.num_completed_sims < numpars.nsims
end

# called after a simulation to update the SimNumberTerminator
@inline function update!(f::SimNumberTerminator, outputter, biopars, numpars)
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

Callback that stops simulating (roughly) when the variance in bound antibodies
becomes sufficiently small.

# Fields
$(FIELDS)
"""
Base.@kwdef mutable struct VarianceTerminator
    """How many simulations have been completed."""
    num_completed_sims::Int
    """The tolerance below which to stop simulating (default = .01)."""
    ssetol::Float64
    """Current estimate of the SSE (default = Inf)."""
    cursse::Float64
end

VarianceTerminator() = VarianceTerminator(0, .01, Inf)

# called before a simulation to see if it should be executed
@inline function isnotdone(vt::VarianceTerminator, biopars, numpars)
    cursse > ssetol
end

# called after a simulation to update the VarianceTerminator
@inline function update!(vt::VarianceTerminator, outputter, biopars, numpars)
    vt.num_completed_sims += 1

    vt.cursse = if vt.num_completed_sims < 15
        Inf
    elseif vt.num_completed_sims > 250
        0.0
    else
        maximum(ouputter.bindcntsse)
    end    

    nothing 
end

# called to reset the VarianceTerminator
@inline function reset!(vt::VarianceTerminator)
    vt.num_completed_sims = 0
    vt.cursse = Inf
    nothing 
end
