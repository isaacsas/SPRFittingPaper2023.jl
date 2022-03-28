@inline function sqdist_periodic(pt1, pt2, L)
    dx = abs(pt1[1] - pt2[1])
    dy = abs(pt1[2] - pt2[2])

    min(dx, L - dx)^2 + min(dy, L - dy)^2 
end

function run_spr_sim(bioparams, numparams)
    @unpack kon,koff,konb,reach,CP,antibodyconcen = bioparams
    @unpack nsims,N,ts_1,tstop,dt,L,initlocs = numparams

    kon  *= antibodyconcen   # convert to per time units
    kc    = 2*koff           # C --> A+B (per time)
    tsave = collect(range(0.0,tstop,step=dt))

    BindCount = zeros(Int,length(tsave))

    reachsq = reach*reach
    rpair = CartesianIndex{2}[]
    Acnt  = zeros(Int,N)
    @inbounds for i = 1:N
        @inbounds for j = (i+1):N
            sqd = sqdist_periodic(view(initlocs,:,i), view(initlocs,:,j), L)
            if sqd <= reachsq
                push!(rpair, CartesianIndex(i,j))
                Acnt[i] += 1
                Acnt[j] += 1
            end
        end
    end

    # now we create a map from a molecule to all possible neighbours
    nhbrs = [Vector{Int}(undef,Acnt[i]) for i in 1:N]
    rxIds = [Vector{Int}(undef,Acnt[i]) for i in 1:N]

    # Create a flipped index version to account for the second ordering of pairs. This will allow us to assume that the first index is 
    # always an A molecule and the second is always a B molecule
    flip(r) = CartesianIndex(r[2],r[1])

    # append the flipped indices, distances, and rates
    append!(rpair, (flip(i) for i in rpair))
    numPairs     = length(rpair)
    halfnumPairs = numPairs ÷ 2

    next  = ones(Int,N)
    ii    = 1
    @inbounds for i=1:numPairs
        ridx             = rpair[i][1]
        idx              = next[ridx]
        nhbrs[ridx][idx] = rpair[i][2]
        rxIds[ridx][idx] = ii
        next[ridx]       = idx+1
        
        ii = (ii == halfnumPairs) ? 1 : (ii + 1)
    end

    # preallocate arrays for current state information
    states      = ones(Int,N,1) # stores the current state of each moleule
    CopyNumbers = zeros(Int,3,length(tsave))
    tvec        = zeros(2*N + numPairs)

     # Now we can run the simulation
    for nr=1:nsims                

        #Initial time
        tc = 0.0

        # Initial particle numbers: All initally time A
        A = N
        B = 0
        C = 0
        states .= 1

        # mean times for each reactions
        τkonb  = 1 / konb
        τkoff = 1 / koff
        τkon  = 1 / kon
        τkc   = 1 / kc
        
        # Saving array
        CopyNumbers .= 0
        CopyNumbers[1,1] = A
        CopyNumbers[2,1] = B
        CopyNumbers[3,1] = C
        sidx             = 2
        tp               = tsave[sidx]

        # Initial reactions times
        tvec .= 0.0
        for i in 1:N
            tvec[i] = τkon*randexp()
        end
        tvec[(N+1):end] .= Inf

        # supposedly this is a faster min heap for storing floats
        times = MutableBinaryHeap{Float64, DataStructures.FasterForward}(tvec)   
        turnoff = false
        twoN = 2*N

        countAB = 0
        countBA = 0
        countABToC = 0
        countCToAB = 0

        @inbounds while tc <= tstop
            if tc >= ts_1 && turnoff == false
                # We turn off the A-->B reaction
                τkon = Inf
                for i=1:N
                    DataStructures.update!(times,i,Inf)
                end
                turnoff = true
            end        

            tnext,rxIdx = top_with_handle(times)

            # update to new time
            tc = tnext

            @inbounds while (tc>tp) && (tp <= tstop)
                CopyNumbers[1,sidx] = A
                CopyNumbers[2,sidx] = B
                CopyNumbers[3,sidx] = C
                tp   += dt
                sidx += 1
            end

            # Update state based on the reaction that has occurred

            if rxIdx <= N
                if tc >= ts_1
                    # then the A->B reaction interacts with the time at which that reaction is turned off
                    molecule         = rxIdx # index for the molecule that just reacted
                    DataStructures.update!(times, molecule, Inf)
                else
                    # A --> B
                    molecule         = rxIdx # index for the molecule that just reacted
                    A                = A-1
                    B                = B+1
                    states[molecule] = 2
                    DataStructures.update!(times, molecule, Inf)
                    DataStructures.update!(times, N + molecule, tc + τkoff * randexp())
                    
                    # Now we check to see if any neighbours can or cannot react this new molecule
                    curNhbrs = nhbrs[molecule] # gives an array of neighbouring molecules
                    curRxs   = rxIds[molecule] # gives the indices of reactions associated to the neighbours
                    
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

            elseif rxIdx <= twoN
                # B --> A
                molecule         = rxIdx - N # index for the molecule that just reacted
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
                curRxs   = rxIds[molecule] # gives the indices of reactions associated to the neighbours
                
                @inbounds for (i,nhbr) in enumerate(curNhbrs)
                    if states[nhbr] == 2 # checks to see that the other molecule is a B
                        # sets a new reaction for the A + B --> C
                        DataStructures.update!(times, twoN + curRxs[i], tc + τkonb * randexp())
                    else
                        # two B's so no reaction can occur
                        DataStructures.update!(times, twoN + curRxs[i], Inf)
                    end
                end


            elseif rxIdx <= (twoN + halfnumPairs)
                # A + B --> C
                curRx = rxIdx - twoN
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
                @inbounds for rxid in rxIds[molecule1]
                    DataStructures.update!(times, twoN + rxid, Inf)
                end
                @inbounds for rxid in rxIds[molecule2]
                    DataStructures.update!(times, twoN + rxid, Inf)
                end                                    
                # turn on possible C --> A + B reaction
                DataStructures.update!(times, twoN + halfnumPairs + curRx, tc + τkc*randexp())

            else
                # C --> A + B
                curRx = rxIdx - (twoN+halfnumPairs)
                A = A + 1
                B = B + 1
                C = C - 1

                DataStructures.update!(times, curRx + twoN + halfnumPairs, Inf)

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
                curRxs   = rxIds[moleculeA]

                @inbounds for (i,nhbr) in enumerate(curNhbrs)
                    if states[nhbr] == 2 # check to see if the neighbour is a B molecule
                        DataStructures.update!(times, twoN + curRxs[i], tc + τkonb * randexp())
                    else
                        DataStructures.update!(times, twoN + curRxs[i] , Inf)
                    end
                end

                curNhbrs = nhbrs[moleculeB]
                curRxs   = rxIds[moleculeB]

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

        @inbounds for i=1:length(tsave)
            BindCount[i] += CopyNumbers[2,i] + CopyNumbers[3,i]
        end
    end

    return BindCount.*(CP/(N*nsims))
end
