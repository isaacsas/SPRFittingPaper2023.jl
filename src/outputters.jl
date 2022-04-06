# structures for saving outputs:

################ SAVING TOTAL BOUND AT EACH TIME ################
struct TotalBoundOutputter
    bindcnt::Vector{Float64}
end

# allows to create
function TotalBoundOutputter(N::Int) 
    TotalBoundOutputter(zeros(N))
end

# this just sums up the amount of B and C at each time
@inline function (o::TotalBoundOutputter)(tsave, copynumbers, biopars, numpars)
    len = size(copynumbers,2)
    for i in 1:len
        o.bindcnt[i] += copynumbers[2,i] + copynumbers[3,i]
    end
    nothing
end

# this is called once all simulations finish to finalize the output
@inline function (o::TotalBoundOutputter)(biopars, numpars)
    @unpack CP = biopars
    @unpack N,nsims = numpars
    o.bindcnt .*= (CP/(N*nsims))
    nothing
end

# this is used to reset the output object prior to a new simulations (i.e. with
# new parameters)
@inline function (o::TotalBoundOutputter)()
    o.bindcnt .= 0
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
