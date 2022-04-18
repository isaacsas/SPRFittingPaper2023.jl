@inline function periodic_dist_sq(pt1, pt2, L) 
    d = 0.0    
    @inbounds for i in eachindex(pt1)
        daxes = abs(pt1[i] - pt2[i])
        d += min(daxes, L - daxes)^2
    end
    d
end

# stores per-simulation specific information
struct NhbrParams{DIM}
   nhbrs::Vector{Vector{Int}} 
   rxids::Vector{Vector{Int}}
   rpair::Vector{CartesianIndex{2}}
   initlocs::Vector{SVector{DIM,Float64}}
   Acnt::Vector{Int}
   next::Vector{Int}
end

# initlocs = DIM by N matrix of the initial particle locations
function NhbrParams(initlocs::Vector{SVector{DIM,Float64}}) where {DIM}
    N     = length(initlocs)
    nhbrs = [Vector{Int}() for _ in 1:N]
    rxids = [Vector{Int}() for _ in 1:N]
    rpair = Vector{CartesianIndex{2}}()
    Acnt  = zeros(Int, N)
    next  = ones(Int, N)
    NhbrParams{DIM}(nhbrs, rxids, rpair, initlocs, Acnt, next)
end

flip(r) = CartesianIndex(r[2],r[1])


# setup molecule initial positions, and then
# for each molecule calculate neighbors within reach
# also setup the reaction index labelling 
function setup_spr_sim!(nhbrpars::NhbrParams{DIM}, biopars, simpars) where {DIM}
    @unpack reach = biopars
    @unpack N,L,resample_initlocs = simpars
    @unpack nhbrs,rxids,rpair,initlocs,Acnt,next = nhbrpars

    if resample_initlocs
        for i = 1:length(initlocs)
            initlocs[i] = rand(SVector{DIM,Float64})
        end
        initlocs .*= L 
    end

    reachsq = reach*reach
    Acnt   .= 0
    empty!(rpair)
    @inbounds for i = 1:N
        @inbounds for j = (i+1):N
            sqd = periodic_dist_sq(initlocs[i], initlocs[j], L)
            if sqd <= reachsq
                push!(rpair, CartesianIndex(i,j))
                Acnt[i] += 1
                Acnt[j] += 1
            end
        end
    end

    # Create a flipped index version to account for the second ordering of pairs. 
    # This will allow us to assume that the first index is 
    # always an A molecule and the second is always a B molecule.
    append!(rpair, (flip(i) for i in rpair))
    numpairs     = length(rpair)
    halfnumpairs = numpairs ÷ 2

    # now we create a map from a molecule to all possible neighbours
    @inbounds for i in 1:N
        (length(nhbrs[i]) != Acnt[i]) && resize!(nhbrs[i], Acnt[i])
        (length(rxids[i]) != Acnt[i]) && resize!(rxids[i], Acnt[i])
    end

    next .= 1
    ii    = 1
    @inbounds for i=1:numpairs
        ridx             = rpair[i][1]
        idx              = next[ridx]
        nhbrs[ridx][idx] = rpair[i][2]
        rxids[ridx][idx] = ii
        next[ridx]       = idx+1
        
        ii = (ii == halfnumpairs) ? 1 : (ii + 1)
    end

    nothing
end

"""
    run_spr_sim!(outputter, biopars::BioPhysParams, simpars::SimParams, 
                 terminator = SimNumberTerminator())

Runs a set of SPR simulations via the particle model for a fixed set of physical
parameters.

Notes:
- Use the passed [`BioPhysParams`](@ref) and [`SimParams`](@ref) to set
  simulation parameters (including number of simulations).
- Use the `outputter` to specify a callback controlling what data is saved.
- Use the `terminator` to specify a callback controlling when the simulation is
  stopped.
"""
function run_spr_sim!(outputter, biopars::BioPhysParams, simpars::SimParams, terminator=SimNumberTerminator())
    @unpack N,tstop_AtoB,tstop,tsave,L,resample_initlocs = simpars

    # don't overwrite the user-provided initlocs
    initlocs = copy(simpars.initlocs)

    # output 
    numsave = length(tsave)

    # setup connectivity info
    nhbrpars = NhbrParams(initlocs)
    setup_spr_sim!(nhbrpars, biopars, simpars)
    @unpack rpair,nhbrs,rxids = nhbrpars
    numpairs     = length(rpair)
    halfnumpairs = numpairs ÷ 2
    twoN         = 2*N
    
    # preallocate arrays for current state information
    states      = ones(Int,N)              # stores the current state of each moleule
    copynumbers = zeros(Int,3,length(tsave))
    tvec        = zeros(2*N + numpairs)
    times       = MutableBinaryHeap{Float64, DataStructures.FasterForward}(tvec)   

     # Now we can run the simulation
    @inbounds while isnotdone(terminator, biopars, simpars)

        # reset simulation specific parameters
        τkon  = 1 / (biopars.kon * biopars.antibodyconcen) # convert to time units
        τkoff = 1 / biopars.koff
        τkonb = 1 / biopars.konb
        τkc   = 1 / (2*biopars.koff)               # C --> A+B (time units)

        # reset initial positions and nhbrs since connectivity changes        
        if resample_initlocs
            setup_spr_sim!(nhbrpars, biopars, simpars)
            numpairs = length(rpair)
            halfnumpairs = numpairs ÷ 2
        end

        #Initial time
        tc = 0.0

        # Initial particle numbers: All initally type A
        A = N; B = 0; C = 0
        states .= 1
        
        # Saving array
        copynumbers     .= 0
        copynumbers[1,1] = A
        copynumbers[2,1] = B
        copynumbers[3,1] = C
        sidx             = 2
        tp               = tsave[sidx]

        # Initial reactions times
        rebuild_times = resample_initlocs && (length(tvec) != (2*N+numpairs))
        rebuild_times && (resize!(tvec, 2*N + numpairs))        
        for i = 1:N
            tvec[i] = τkon*randexp()
        end
        tvec[(N+1):end] .= Inf

        if rebuild_times
            times = MutableBinaryHeap{Float64, DataStructures.FasterForward}(tvec)   
        else 
            for (i,tval) in pairs(tvec)
                DataStructures.update!(times, i, tval)
            end
        end
        turnoff = false

        @inbounds while tc <= tstop
            if tc >= tstop_AtoB && turnoff == false
                # We turn off the A-->B reaction
                τkon = Inf
                for i=1:N
                    DataStructures.update!(times,i,Inf)
                end
                turnoff = true
            end        

            tnext,rxidx = top_with_handle(times)

            # update to new time
            tc = tnext

            @inbounds while (tc>tp) && (sidx <= numsave)
                copynumbers[1,sidx] = A
                copynumbers[2,sidx] = B
                copynumbers[3,sidx] = C                
                sidx += 1
                tp = (sidx <= numsave) ? tsave[sidx] : (tstop+eps(tstop))
            end

            # Update state based on the reaction that has occurred

            if rxidx <= N
                if tc >= tstop_AtoB
                    # then the A->B reaction interacts with the time at which that reaction is turned off
                    molecule         = rxidx # index for the molecule that just reacted
                    DataStructures.update!(times, molecule, Inf)
                else
                    # A --> B
                    molecule         = rxidx # index for the molecule that just reacted
                    A                = A-1
                    B                = B+1
                    states[molecule] = 2
                    DataStructures.update!(times, molecule, Inf)
                    DataStructures.update!(times, N + molecule, tc + τkoff * randexp())
                    
                    # Now we check to see if any neighbours can or cannot react this new molecule
                    curNhbrs = nhbrs[molecule] # gives an array of neighbouring molecules
                    curRxs   = rxids[molecule] # gives the indices of reactions associated to the neighbours
                    
                    @inbounds for (i,nhbr) in enumerate(curNhbrs)
                        if states[curNhbrs[i]] == 1 # checks to see that the other molecule is an A
                            # sets a new reaction for the A + B --> C
                            DataStructures.update!(times, twoN + curRxs[i], tc + τkonb*randexp())
                        else
                            # two B's so no reaction can occur
                            DataStructures.update!(times, twoN + curRxs[i],  Inf)
                        end
                    end

                end

            elseif rxidx <= twoN
                # B --> A
                molecule         = rxidx - N # index for the molecule that just reacted
                A                = A+1
                B                = B-1
                states[molecule] = 1
                DataStructures.update!(times, N + molecule, Inf)
                if turnoff==false
                    DataStructures.update!(times, molecule, tc + τkon * randexp())
                else
                    DataStructures.update!(times, molecule, Inf)
                end

                # Now we check to see if any neighbours can or cannot react this new molecule
                curNhbrs = nhbrs[molecule] # gives an array of neighbouring molecules
                curRxs   = rxids[molecule] # gives the indices of reactions associated to the neighbours
                
                @inbounds for (i,nhbr) in enumerate(curNhbrs)
                    if states[nhbr] == 2 # checks to see that the other molecule is a B
                        # sets a new reaction for the A + B --> C
                        DataStructures.update!(times, twoN + curRxs[i], tc + τkonb * randexp())
                    else
                        # two B's so no reaction can occur
                        DataStructures.update!(times, twoN + curRxs[i], Inf)
                    end
                end


            elseif rxidx <= (twoN + halfnumpairs)
                # A + B --> C
                curRx = rxidx - twoN
                A = A-1
                B = B-1
                C = C+1

                # AB configuration
                molecule1 = rpair[curRx][1]
                molecule2 = rpair[curRx][2]

                # swap if BA configuration
                if states[molecule1] == 2
                    molecule1,molecule2 = molecule2,molecule1
                end

                states[molecule1] = 0
                states[molecule2] = 0
                
                # turn off A <--> B reactions 
                DataStructures.update!(times, molecule1, Inf)
                DataStructures.update!(times, molecule2, Inf)
                DataStructures.update!(times, molecule1+N, Inf)
                DataStructures.update!(times, molecule2+N, Inf)
                            
                # turn off all A+B --> C reactions
                @inbounds for rxid in rxids[molecule1]
                    DataStructures.update!(times, twoN + rxid, Inf)
                end
                @inbounds for rxid in rxids[molecule2]
                    DataStructures.update!(times, twoN + rxid, Inf)
                end                                    
                # turn on possible C --> A + B reaction
                DataStructures.update!(times, twoN + halfnumpairs + curRx, tc + τkc*randexp())

            else
                # C --> A + B
                curRx = rxidx - (twoN+halfnumpairs)
                A = A + 1
                B = B + 1
                C = C - 1

                DataStructures.update!(times, curRx + twoN + halfnumpairs, Inf)

                pAB = rand() # random number to decide which molecule to place first, i.e. which tether is freed
                moleculeA = rpair[curRx][1]
                moleculeB = rpair[curRx][2]

                # swap the molecules in this case
                if pAB >= 0.5
                    moleculeA,moleculeB = moleculeB,moleculeA
                end

                # Then A molecule is placed first
                states[moleculeA] = 1
                states[moleculeB] = 2   

                if turnoff==false
                    DataStructures.update!(times, moleculeA, tc + τkon * randexp())
                else
                    DataStructures.update!(times, moleculeA, Inf)
                end
                DataStructures.update!(times, N + moleculeA, Inf)
                DataStructures.update!(times, moleculeB, Inf)
                DataStructures.update!(times, N + moleculeB, tc + τkoff * randexp())
                
                curNhbrs = nhbrs[moleculeA]
                curRxs   = rxids[moleculeA]

                @inbounds for (i,nhbr) in enumerate(curNhbrs)
                    if states[nhbr] == 2 # check to see if the neighbour is a B molecule
                        DataStructures.update!(times, twoN + curRxs[i], tc + τkonb * randexp())
                    else
                        DataStructures.update!(times, twoN + curRxs[i] , Inf)
                    end
                end

                curNhbrs = nhbrs[moleculeB]
                curRxs   = rxids[moleculeB]

                @inbounds for (i,nhbr) in enumerate(curNhbrs)
                    # check to see if the neighbour is an A molecule
                    if (states[nhbr] == 1) && (nhbr != moleculeA)
                        DataStructures.update!(times, twoN + curRxs[i], tc + τkonb * randexp())                    
                    elseif states[nhbr] != 1
                        DataStructures.update!(times, twoN + curRxs[i], Inf)
                    end
                end                 
            end
        end

        # save at remaining output times (which must be greater than stopping
        # time) this assumes solution is constant past tstop!
        @inbounds while sidx <= numsave
            copynumbers[1,sidx] = A
            copynumbers[2,sidx] = B
            copynumbers[3,sidx] = C
            sidx += 1
        end

        # save per simulation output
        outputter(tsave, copynumbers, biopars, simpars)

        # update terminator callback
        update!(terminator, outputter, biopars, simpars)
    end

    # post simulation processing of output
    outputter(biopars, simpars)
    nothing
end
