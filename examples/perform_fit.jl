using SPRFitting
using LinearAlgebra, Interpolations, LossFunctions, BlackBoxOptim
using DataFrames, CSV, DelimitedFiles, XLSX
using Plots, JLD


############ INPUT ############ 
experiment_name = "Test_110322"
# Note this vector may need to be changed for different experiments as the antigen concentrations are changed
#global ABC = [18.75, 75., 300.] #nM
# edit: the code takes ABC from header in input file. 
LUTname = "HighKon" # choose "LowKon" or "HighKon"

save_curves = true
visualise = true

############ FUNCTIONS ############

function InterpolateFunctionError(param)
    # This function takes a set of parameters and interpolates a simulated kinetics curve from the Look Up Table 
    # that we made in Tutorial 1, and then calculate the error between the interpolated curve and the experimental data
    # using a sum of squares.
    #global itp
    #global n_gp
    #global RD

    n_gp = size(RD)[2]
    
    # Declare the lookuptable bounds as global
    #global p1_min_LUT,p1_range_LUT
    #global p2_min_LUT,p2_range_LUT
    #global p3_min_LUT,p3_range_LUT
    #global Lmin_LUT,L_range_LUT


    # Rescale the parameters into their indexings for the interpolant.
    # Anna: What is the 29 in this code??
    q2 = (param[2]-p2_min_LUT)*(29/(p2_range_LUT))+1
    q3 = (param[3]-p3_min_LUT)*(29/(p3_range_LUT))+1
    qL = (param[4]-Lmin_LUT)*(29/(L_range_LUT))+1

    e1 = zeros(length(J),n_gp)
    #ABC = [0.2,0.39,0.78,1.56,3.1,6.3,12.5,25] # antibody concentrations in nM
    for (jx,j) in enumerate(J)
            
        k1 = param[1] + log10(ABC[j]/ABC[1]) # rescale the on rate to account for changing concentraions of antibody

        q1 = (k1-p1_min_LUT)*(29/(p1_range_LUT))+1
        #qL = (param[j+8]-5)*(14/30)+1
        
        for i=1:n_gp
            if isnan.(Times[i]) == 0
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

function SimulateFunction(param)

    # Set up the parameters that wont change in each simulation
    AT = 1.0 # antibody concentration
    #biophysical on and off rates
    #k1 = 0.079
    #km1 = 0.174
    #k2 = 0.05
    # D and conversion factor
    convfactor = 6.023 * 1e-7

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
    tfinal = 600 # make less rigid?
    ts_1 = 150 # time at which we shut off the A --> B reaction
    dt = 1 # time at which we save the data
    t = 1:1:(tfinal+1) # time vec for interpolation
    initLocs = rand(3,N)*L   # moved this here...

    Simulate_AntibodyModel(kc,100,kon,koff,lam,LAB,N,L,ts_1,tfinal,dt,tsave,CP)
end

function createData(file)
    # reads Aligned Data and returns Time, Data and Antibody concentration
    Data, Header = readdlm(joinpath(RAWDIR, file),',', header=true)
    Abconc = parse.(Float64,Header[2:end])

    Times = transpose(Data[:,1])
    RefData = transpose(Data[:, 2:end])
    println(size(Data), size(RefData))
    #Ignore = transpose(GapsData)

    RD = ifelse.(isnan.(RefData[1:end,:]), 0,RefData[1:end,:])

    return Times, RefData, RD, Abconc
end

function getLUT(LUTname)
    # sets parameters for Lookup table used
    # Clear Garbage
    LUT_Data=0
    itp=0
    GC.gc()

    # choose which LUT (currently two)
    if LUTname == "LowKon"
        LUT_Data = load(joinpath(BASEDIR,"Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld"))["FirstMoment"]
    elseif LUTname == "HighKon"
        LUT_Data = load(joinpath(BASEDIR,"Surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld"))["FirstMoment"]
    end

    global itp = interpolate(LUT_Data,BSpline(Linear()))

    # Release the LUT from memory
    LUT_Data = nothing
    GC.gc()

    # bounds on the reach
    global Lmin_LUT = 2
    global Lmax_LUT = 35
    global L_range_LUT = Lmax_LUT-Lmin_LUT


    # choose which LUT (currently two)
    if LUTname == "LowKon"
        # this is the bounds of the LookUpTable used here and the table 
        # 'CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27'
        # which is found in Codes/Surrogates
        global p1_min_LUT = -5.0
        global p1_max_LUT = -0.0 
        global p1_range_LUT = p1_max_LUT - p1_min_LUT

    elseif LUTname == "HighKon"
        # this is the bounds of the LookUpTable
        #' CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Jan27'
        # which is found in Codes/Surrogates

        global p1_min_LUT = -3.0
        global p1_max_LUT = 2.0 
        global p1_range_LUT = p1_max_LUT - p1_min_LUT

    end
    
    global p2_min_LUT = -4.0
    global p2_max_LUT = -1.0 
    global p2_range_LUT = p2_max_LUT - p2_min_LUT

    global p3_min_LUT = -3.0
    global p3_max_LUT = 1.5 
    global p3_range_LUT = p3_max_LUT - p3_min_LUT
    

    # Prior bounds for the search, these might be more restrictive than the bounds upon which the LUT was made
    # in order for a more efficient localised search, and then can be widened as necessary to check there are no
    # better global optima

    Lmin = 2.0 # LB for molecular reach
    Lmax = 35.0 # UC for molecular reach
    p1_min = -5.0 # LB for log10(k_on)
    p1_max = -1.25
   # p1_max = -2.5 # UB for log10(k_on) NB: its lower than the bounds on the LUT as we need to be able to search all the antibody concentrations at once.
    p2_min = -4.0 # LB for log10(k_off)
    p2_max = -1.0 # UB for log10(k_off)
    p3_min = -3.0 # LB for log(k_on,b)
    p3_max = 1.0  # UB for log10(k_on,b)

    CP_min = log10(10) # LB for the constant of proportionality
    CP_max = log10(100000) # UB for the constant of proportionality

    custombounds = [(p1_min,p1_max),(p2_min,p2_max),(p3_min,p3_max),(Lmin,Lmax),(CP_min,CP_max)]
    return custombounds
end

function visualiseFit(result, Times, RefData,filename)
    # plots data and simulated curves and saves it in OUTDIR
    Params=copy(result.method_output.population[1].params)
    # @show "Is this code being run?"
    #n_dil = size(RefData)[1]
    n_dil = size(J)[1]
    fig1 = plot(Times',RefData[1,:],label="",color="black")
    for j= 2:n_dil
        plot!(Times',RefData[j,:],label="",color="black")
        # if j ==1
        #     fig1 = plot(Times',RefData[j,:],label="",color="black")
        # else
        #     plot!(Times',RefData[j,:],label="",color="black")
        # end
    end
    

    timepoints=0:600#Times[end]
    for j= 1:n_dil
        params = [Params[1]+log10(ABC[j]/ABC[1]),Params[2],Params[3],Params[4],Params[5]]
        Y=SimulateFunction(params)
        plot!(timepoints,Y[1,:],label="")
    end

    xlabel!("time")
    ylabel!("RU")
    savefig(fig1,  joinpath(OUTDIR, filename * "_Fitted.png"))
end

function saveFit(result, Times, RefData, filename)
    # saves data and simulated data in OUTDIR
    Params=copy(result.method_output.population[1].params)
    n_dil = size(J)[1] #size(RefData)[1]
    n_gp = size(RefData)[2]
    timepoints=0:600 #Times[end]
    SaveData = zeros(Float64,length(timepoints),2*n_dil+1)
    SaveData[:,1] .= timepoints
    

    for j= 1:n_dil

        # Assign the aligned data to match timepoints
        nstart = Int(Times[1])
        SaveData[1:nstart,2*j] .= NaN # Times starts at 2 -5 sec. 
        for n=1:n_gp # length(timepoints) - nstart
            SaveData[nstart+n,2*j] = RefData[j,n] # already contains Nan
        end

        # Simulate the model and include the averaged kinetics data
        params = [Params[1]+log10(ABC[j]/ABC[1]),Params[2],Params[3],Params[4],Params[5]]
        Y=SimulateFunction(params) # this function assumes time = 600 sec
        SaveData[:,2*j+1] .= Y[1,:]
    end

    # Create a DataFrame with the columns
    df = DataFrame(SaveData,:auto)
    EH = [Symbol("Experimental Data $(ABC[j]) nM") for j in 1:n_dil]
    MH = [Symbol("Model Data $(ABC[j]) nM") for j in 1:n_dil]
    CH = [Z[i] for i=1:n_dil for Z in [EH,MH]]

    Headers = [Symbol("times")]
    for h in CH
        push!(Headers,h)
    end
    
    # Rename the headers of each columns
    rename!(df,Headers)

    # Write the DataFrame to an xlsx file
    XLSX.writetable(joinpath(OUTDIR, filename*"_Fitted.xlsx"),collect(eachcol(df)),names(df), overwrite=true)

    ##### save parameters
    bs = best_candidate(result)
    bf = best_fitness(result)
    open(joinpath(OUTDIR, filename*"_Fitted.txt"), "w+") do file
        println(file, filename, "\n")
        println(file, "Best candidate found (kon, koff, konb, reach, CP): ", bs)
        println(file, "Fitness: ", bf)
    end
end


############ RUN SCRIPT ############

# Directories 
BASEDIR = joinpath(@__DIR__, "figures and data")
@show BASEDIR
RAWDIR = joinpath(BASEDIR, "Experiments","Aligned", experiment_name)
OUTDIR = joinpath(BASEDIR, "Experiments", "Fitted", experiment_name)

# this should create the directory if it doesn't exist
mkpath(OUTDIR)

# Loop through files and do the fitting
allfiles = readdir(RAWDIR)
custombounds = getLUT(LUTname)

# for file in allfiles
#     if occursin(r"^Data_", file) == true
#         filename = replace(replace(file,r"^Data_" => "" ), r".csv$" => "")
#         println(filename, file)

#         Times, RefData, RD, ABC = createData(file)
#         print(ABC)
#         J =[1:1:size(RD)[1];]
        
#         # I need to set these variables global, because I don't know how to feed them into the bboptimize function
#         # not the nicest solution but I don"t know how to do it otherwise
#         global RD
#         global J
#         global Times
#         global ABC
#         result = bboptimize(InterpolateFunctionError; SearchRange=custombounds, NumDimensions=5, Method=:xnes,MaxSteps=5000, TraceMode=:compact, TraceInterval=10.0)#,Tracer=:silent)
#         if visualise == true
#             visualiseFit(result, Times, RefData, filename)
#         end
#         if save_curves == true
#             saveFit(result, Times, RefData, filename)
#         end
#     end
# end


