using BlackBoxOptim, Random, Distributions, StatsBase, DataFrames, XLSX
using DataStructures, LinearAlgebra, Plots, JLD, Interpolations, LossFunctions


# This function takes a set of parameters and interpolates a simulated kinetics
# curve from the Look Up Table that we made in Tutorial 1, and then calculate
# the error between the interpolated curve and the experimental data using a sum
# of squares.
function InterpolateFunctionError(param)

    global itp
    global n_gp
    global RD
    
    # Declare the lookuptable bounds as global
    global p1_min_LUT,p1_range_LUT
    global p2_min_LUT,p2_range_LUT
    global p3_min_LUT,p3_range_LUT
    global Lmin_LUT,L_range_LUT


    # Rescale the parameters into their indexings for the interpolant.
    q2 = (param[2]-p2_min_LUT)*(29/(p2_range_LUT))+1
    q3 = (param[3]-p3_min_LUT)*(29/(p3_range_LUT))+1
    qL = (param[4]-Lmin_LUT)*(29/(L_range_LUT))+1

        e1 = zeros(length(J),n_gp)
        # Note this vector may need to be changed for different experiments as the antigen concentrations are changed
        ABC = [0.2,0.39,0.78,1.56,3.1,6.3,12.5,25] # antibody concentrations in nM
    for (jx,j) in enumerate(J)
            
            k1 = param[1] + log10(ABC[j]/ABC[1]) # rescale the on rate to account for changing concentraions of antibody

            q1 = (k1-p1_min_LUT)*(29/(p1_range_LUT))+1
            #qL = (param[j+8]-5)*(14/30)+1
            
            for i=1:n_gp
                if Ignore[j,i]==0
                    e1[jx,i] = itp(q1,q2,q3,qL,Int(Times[i]))*10^(param[5]) 
                    #RD[j,i] = RefData[j,i]
                    #Weights[j,i] = 1 ./ RefData[j,i]^2
                end  
            end
        end
        # @show size(e1),size(RD)
    error = sqrt(value(L2DistLoss(),e1[:],RD[:],AggMode.Sum()))
   
    return error
end


@inline function sqdist_periodic(pt1, pt2, L)
    dx = abs(pt1[1] - pt2[1])
    dy = abs(pt1[2] - pt2[2])
    dz = abs(pt1[3] - pt2[3])

    min(dx, L - dx)^2 + min(dy, L - dy)^2 + min(dz, L - dz)^2
end


function SimulateFunction(param)

    # Set up the parameters that wont change in each simulation
    AT = 1.0 # antibody concentration
    #biophysical on and off rates
    #k1 = 0.079
    #km1 = 0.174
    #k2 = 0.05
    # D and conversion factor
    D = 15
        convfactor = 6.023 * 1e-7
    # the molecular concentration of the second-step reaction
    sigD = 3/(4*pi*D^3*convfactor)
    # molecular reach of the second-step reaction
    LAB = param[4]
    # constant of proportionality
    CP = 10^(param[5])

    # Convert biophysical parameters into model parameters
    kon = 10^(param[1]) * AT       # A --> B
    koff = 10^(param[2])       # B --> A
    lam =  10^(param[3]) # A+B --> C    --- Makes it orthogonal as now param[1] = k2 / km2
    kc = 10^(param[2])*2

    # Simulation Control
    ST_uM = 500/149/26795*1000000 # Maps 500 RU to a concentration in uM
    ST = ST_uM * convfactor # converts from uM to molecules per nm^3
    N = 1000 # number of particles (tethers)
    L = (N/ST)^(1/3) # so we choose both the concentration and the number of particles and find the appropriate size domain
    tfinal = 600
    ts_1 = 150 # time at which we shut off the A --> B reaction
    dt = 1 # time at which we save the data
    t = 1:1:(tfinal+1) # time vec for interpolation

    function Simulate_AntibodyModel(kc,nsim,kon,koff,lam,LAB,N,L,ts_1,tfinal,dt,tsave,CP)

        BindCount = zeros(Int,1,length(tsave))
        for nr=1:nsim
            
            initLocs = rand(3,N)*L

            LABsq = LAB*LAB
            rpair = CartesianIndex{2}[]
            Acnt  = zeros(Int,N)
            @inbounds for i = 1:N
                @inbounds for j = (i+1):N
                    sqd = sqdist_periodic(view(initLocs,:,i), view(initLocs,:,j), L)
                    if sqd <= LABsq
                        push!(rpair, CartesianIndex(i,j))
                        Acnt[i] += 1
                        Acnt[j] += 1
                    end
                end
            end
            # Create a flipped index version to account for the second ordering of pairs. This will allow us to assume that the first index is 
            # always an A molecule and the second is always a B molecule
            flip(r) = CartesianIndex(r[2],r[1])

            # append the flipped indices, distances, and rates
            append!(rpair, (flip(i) for i in rpair))
            numPairs     = length(rpair)
            halfnumPairs = numPairs ÷ 2

            # now we create a map from a molecule to all possible neighbours
            nhbrs = [Vector{Int}(undef,Acnt[i]) for i in 1:N]
            rxIds = [Vector{Int}(undef,Acnt[i]) for i in 1:N]

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

            # Now we can run the simulation

            #Initial time
            tc = 0.0

            # Initial particle numbers: All initally time A
            A = N
            B = 0
            C = 0
            states = ones(Int,N,1) # stores the current state of each moleule

            # mean times for each reactions
            τlam  = 1 / lam
            τkoff = 1 / koff
            τkon  = 1 / kon
            τkc   = 1 / kc
            
            # Saving array
            CopyNumbers = zeros(Int,3,length(tsave))
            CopyNumbers[1,1] = A
            CopyNumbers[2,1] = B
            CopyNumbers[3,1] = C
            sidx             = 2
            tp               = tsave[sidx]

            # Initial reactions times
            tvec = zeros(2*N + numPairs)
            foreach(i -> tvec[i] = τkon*randexp(), 1:N)
            tvec[(N+1):end] .= Inf

            # supposedly this is a faster min heap for storing floats
            times = MutableBinaryHeap{Float64, DataStructures.FasterForward}(tvec)   
            turnoff = false
            twoN = 2*N

            countAB = 0
            countBA = 0
            countABToC = 0
            countCToAB = 0

            @inbounds while tc <= tfinal


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

                @inbounds while (tc>tp) && (tp <= tfinal)
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
                                DataStructures.update!(times, twoN + curRxs[i], tc + τlam*randexp())
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
                            DataStructures.update!(times, twoN + curRxs[i], tc + τlam * randexp())
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
                    DataStructures. update!(times, N + moleculeA, Inf)
                    DataStructures.update!(times, moleculeB, Inf)
                    DataStructures.update!(times, N + moleculeB, tc + τkoff * randexp())
                    
                    curNhbrs = nhbrs[moleculeA]
                    curRxs   = rxIds[moleculeA]

                    @inbounds for (i,nhbr) in enumerate(curNhbrs)
                        if states[nhbr] == 2 # check to see if the neighbour is a B molecule
                            DataStructures.update!(times, twoN + curRxs[i], tc + τlam * randexp())
                        else
                            DataStructures.update!(times, twoN + curRxs[i] , Inf)
                        end
                    end

                    curNhbrs = nhbrs[moleculeB]
                    curRxs   = rxIds[moleculeB]

                    @inbounds for (i,nhbr) in enumerate(curNhbrs)
                        # check to see if the neighbour is an A molecule
                        if (states[nhbr] == 1) && (nhbr != moleculeA)
                            DataStructures.update!(times, twoN + curRxs[i], tc + τlam * randexp())                    
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
            BindCount.*(CP/(N*nsim))
    end
    tsave = collect(0:dt:tfinal)
    Simulate_AntibodyModel(kc,100,kon,koff,lam,LAB,N,L,ts_1,tfinal,dt,tsave,CP)
end

J = [1,2,3,4,5,6]
#RD = makeRD(RefData,Ignore)

#SPRData = readdlm("AlignedSPRData/SPR_07_02_2022_FD5D_17uM.csv",',')
#GapsData = readdlm("AlignedSPRData/Gaps_07_02_2022_FD5D_17uM.csv",',')
SPRData = readdlm("/Users/omer/Dropbox/Collaborations/2022 - Bivalent Antibody Project/21-02-2022-Dans_Tutorials_and_Codes/Codes/Parameter Fitting/AlignedSPRData/SPR_07_02_2022_FD5D_17uM.csv",',')
GapsData = readdlm("/Users/omer/Dropbox/Collaborations/2022 - Bivalent Antibody Project/21-02-2022-Dans_Tutorials_and_Codes/Codes/Parameter Fitting/AlignedSPRData/Gaps_07_02_2022_FD5D_17uM.csv",',')



#CurveID = 3 # which curve in the dataset we fit to

n_gp = size(SPRData)[1] # number of grid points in the SPR data
#n_gp=400
RefData = zeros(Float64,7,n_gp)
Times = zeros(Float64,1,n_gp)
Ignore = zeros(Int64,7,n_gp)

for i=1:n_gp
    Times[i] = SPRData[i,1]
end

for j=1:7
    for i=1:n_gp
        RefData[j,i] = SPRData[i,j+1]
        Ignore[j,i] = GapsData[i,j]
    end
end
RD = makeRD(RefData,Ignore)
# Clear Garbage
LUT_Data=0
itp=0
GC.gc()

#RefData = SimulateFunction([-0.12,-0.34,-0.73,-2.5,14.5,0.0])

# True values p1 = 0.0, p2 = -0.30103, p3 = -0.69897, p4 = -2.60206, p5 = 15, p6 = 0

# LUT_Data = load("/home/daniel/Dropbox (BOSTON UNIVERSITY)/Antibody_Project/21-02-2022-Dans_Tutorials_and_Codes/Codes/Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld")["FirstMoment"]
# LUT_Data = load("/home/daniel/Documents/Antibody Project/LookUpTables/NewTablesCovid/Separated_LUTs/CombinedSurrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld")["FirstMoment"]
LUT_Data = load("/Users/omer/Dropbox/Collaborations/2022 - Bivalent Antibody Project/21-02-2022-Dans_Tutorials_and_Codes/Codes/Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld")["FirstMoment"]

itp = interpolate(LUT_Data,BSpline(Linear()))
# Release the LUT from memory
LUT_Data=0
LUT_Data_Old=0
GC.gc()

# bounds on the reach
Lmin_LUT = 2
Lmax_LUT = 35
L_range_LUT = Lmax_LUT-Lmin_LUT

# this is the bounds of the LookUpTable used here and the table 
# 'CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27'
# which is found in Codes/Surrogates
p1_min_LUT = -5.0
p1_max_LUT = -0.0 
p1_range_LUT = p1_max_LUT - p1_min_LUT

# this is the bounds of the LookUpTable
#' CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Jan27'
# which is found in Codes/Surrogates

# p1_min_LUT = -3.0
# p1_max_LUT = 2.0 
# p1_range_LUT = p1_max_LUT - p1_min_LUT

p2_min_LUT = -4.0
p2_max_LUT = -1.0 
p2_range_LUT = p2_max_LUT - p2_min_LUT

p3_min_LUT = -3.0
p3_max_LUT = 1.5 
p3_range_LUT = p3_max_LUT - p3_min_LUT
 

# Prior bounds for the search, these might be more restrictive than the bounds upon which the LUT was made
# in order for a more efficient localised search, and then can be widened as necessary to check there are no
# better global optima

Lmin = 2.0 # LB for molecular reach
Lmax = 35.0 # UC for molecular reach
p1_min = -5.0 # LB for log10(k_on)
p1_max = -2.5 # UB for log10(k_on) NB: its lower than the bounds on the LUT as we need to be able to search all the antibody concentrations at once.
p2_min = -4.0 # LB for log10(k_off)
p2_max = -1.0 # UB for log10(k_off)
p3_min = -3.0 # LB for log(k_on,b)
p3_max = 1.0  # UB for log10(k_on,b)

CP_min = log10(10) # LB for the constant of proportionality
CP_max = log10(100000) # UB for the constant of proportionality

custombounds = [(p1_min,p1_max),(p2_min,p2_max),(p3_min,p3_max),(Lmin,Lmax),(CP_min,CP_max)]
result = bboptimize(InterpolateFunctionError; SearchRange=custombounds, NumDimensions=5, Method=:xnes,MaxSteps=5000,Tracer=:verbose)

visualise = true

save_curves = false

if visualise == true
    Params=copy(result.method_output.population[1].params)
    # @show "Is this code being run?"

    for (jx,j) in enumerate(J)
        if jx==1
            display(plot(Times',RefData[j,:],label="",color="black"))
        else
            display(plot!(Times',RefData[j,:],label="",color="black"))
        end
    end
    

    # AGC = [3.90625,7.8125,15.625,31.25,62.5,125,250,500,1000]   # antibody concentrations in nM
    # AGC = [2.34,4.69,9.38,18.75,37.5,75,150,300]
    AGC = [9.38,18.75,37.5,75,150,300]

    timepoints=0:600
    for j in J
        params = [Params[1]+log10(AGC[j]/AGC[1]),Params[2],Params[3],Params[4],Params[5]]
        Y=SimulateFunction(params)
        display(plot!(timepoints,Y[1,:],label=""))
    end

    xlabel!("time")
    ylabel!("RU")
end

if save_curves == true
    Params=copy(result.method_output.population[1].params)
    
    timepoints=0:600
    SaveData = zeros(Float64,length(timepoints),2*length(J)+1)
    SaveData[:,1] .= timepoints
    
    AGC = [9.38,18.75,37.5,75,150,300]

    for (jx,j) in enumerate(J)

        # Assign the aligned data to match timepoints
        nstart = Int(Times[1])
        SaveData[1:nstart,2*jx] .= NaN
        for n=1:n_gp
            if Ignore[j,n]==0
                SaveData[nstart+n,2*jx] = RefData[j,n]
            else
                SaveData[nstart+n,2*jx] = NaN
            end
        end

        # Simulate the model and include the averaged kinetics data
        params = [Params[1]+log10(AGC[j]/AGC[1]),Params[2],Params[3],Params[4],Params[5]]
        Y=SimulateFunction(params)
        SaveData[:,2*jx+1] .= Y[1,:]
    end

    # Create a DataFrame with the columns
    df = DataFrame(SaveData,:auto)
    EH = [Symbol("Experimental Data $(AGC[j]) nM") for j in J]
    MH = [Symbol("Model Data $(AGC[j]) nM") for j in J]
    CH = [Z[i] for i=1:length(J) for Z in [EH,MH]]

    Headers = [Symbol("times")]
    for h in CH
        push!(Headers,h)
    end
    
    # Rename the headers of each columns
    rename!(df,Headers)

    # Write the DataFrame to an xlsx file
    XLSX.writetable("CurveData/07-02-2022-FD5D-AGC-17uM.xlsx",collect(eachcol(df)),names(df))
end




