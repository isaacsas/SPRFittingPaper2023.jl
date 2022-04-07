using SPRFitting
using LinearAlgebra
using Plots

# to convert 500 RU supposedly to μM
const SIMULATION_ANTIGEN_CONCEN = 6.023*.1*500/149/26795

############ INPUT ############ 
experiment_name = "Test_110322"

# file directories
#BASEDIR = joinpath(@__DIR__, "figures and data")
BASEDIR = "/Users/isaacsas/Desktop/anna_fitting/BivalentAntibody"
RAWDIR = joinpath(BASEDIR, "Experiments","Aligned", experiment_name)
OUTDIR = joinpath(BASEDIR, "Experiments", "Fitted_Sam", experiment_name)

# high kon
lutfile = "/Users/isaacsas/Dropbox/Collaborations/Omer Dushek/2022 - Bivalent Antibody Project/21-02-2022-Dans_Tutorials_and_Codes/Codes/Surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld"
#lutfile = joinpath(BASEDIR,"Surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld")
logkon_range = (-3.0,2.0)    

# low kon
#lutfile  = joinpath(BASEDIR,"Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld")
#logkon_range = (-5.0,-0.0)

# output control
save_curves = true
visualise   = true

# surrogate parameter ranges (log space except reach)
logkoff_range = (-4.0, -1.0)
logkonb_range = (-3.0, 1.5)
reach_range   = (2.0, 35.0)
logCP_range   = (1.0, 5.0)

# optimizer parameter ranges (log space except reach)
#logkon_optrange  = (-5.0, -1.25)  # or -2.5 in old file
logkon_optrange = (-3.0, -1.25)
logkoff_optrange = logkoff_range
logkonb_optrange = (-3.0, 1.0)
reach_optrange   = reach_range
logCP_optrange   = logCP_range

# forward simulation parameters
simpars = SimParams(SIMULATION_ANTIGEN_CONCEN; 
                    N = 1000,           # number of particles
                    tstop = 600.0,      # time to end simulations
                    tstop_AtoB = 150.0, # time to shut off A --> B
                    dt = 1.0,           # time frequency to save 
                    DIM=3,              # use a cubic domain
                    nsims = 100)        # number of simulations to run

############################ END INPUT ###############################

#################### INTERNAL PARAMETERS #############################

# the surrogate's parameter ranges (all log space except reach)
logpar_ranges = SurrogateRanges(; logkon_range, logkoff_range, logkonb_range, reach_range, logCP_range)

# check optimizer ranges are subset of surrogate ranges
checkrange(rsur,ropt) = (rsur[1] <= ropt[1] <= ropt[2] <= rsur[2])
checkrange(logkon_range, logkon_optrange) || error("logkon_range not subset of logkon_optrange")
checkrange(logkoff_range, logkoff_optrange) || error("logkoff_range not subset of logkoff_optrange")
checkrange(logkonb_range, logkonb_optrange) || error("logkonb_range not subset of logkonb_optrange")
checkrange(reach_range, reach_optrange) || error("reach_range not subset of reach_optrange")
checkrange(logCP_range, logCP_range) || error("logCP_range not subset of logCP_optrange")

# the optimizer's parameter ranges (all log space except reach)
optpar_ranges = [logkon_optrange,logkoff_optrange,logkonb_optrange,reach_optrange,logCP_optrange]

# build the surrogate 
surrogate = Surrogate(logpar_ranges, lutfile)

################## END INTERNAL PARAMETERS ############################

############ RUN SCRIPT ############

# this should create the directory if it doesn't exist
mkpath(OUTDIR)

# Loop through files and do the fitting
allfiles = readdir(RAWDIR)

for file in allfiles
    if occursin(r"^Data_", file) == true
        filename = replace(replace(file,r"^Data_" => "" ), r".csv$" => "")

        # get the alinged SPR data
        aligneddat = get_aligned_data(joinpath(RAWDIR, file))

        # find the best fit parameters
        bbopt_output = fit_spr_data(surrogate, aligneddat, optpar_ranges)

        if visualise
            figfile = joinpath(OUTDIR, filename * "_fit_curves.png")
            visualisefit(bbopt_output, aligneddat, simpars, figfile)
        end
        if save_curves
            curvefile = joinpath(OUTDIR, filename)
            savefit(bbopt_output, aligneddat, simpars, curvefile)
        end
    end
end


