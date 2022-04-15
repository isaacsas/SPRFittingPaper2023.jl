# structures for saving outputs:

################ SAVING TOTAL BOUND AT EACH TIME ################
struct TotalBoundOutputter
    """Average amount bound"""
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

# this will be used by the simulator to store output as it goes
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

# structures for determining if need more simulations or not

# stops simulations when the desired number of simulations is reached
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
