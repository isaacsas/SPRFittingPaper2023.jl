using SPRFitting
using LinearAlgebra
using Plots

############ INPUT ############

experiment_name = "Test_110322"

# file directories
#BASEDIR = joinpath(@__DIR__, "figures and data")
BASEDIR = "/Users/isaacsas/Desktop/anna_fitting/BivalentAntibody"
RAWDIR = joinpath(BASEDIR, "Experiments","Aligned", experiment_name)
OUTDIR = joinpath(BASEDIR, "Experiments", "Fitted_Sam", experiment_name)

# SPR antigen concentration in μM
antigenconcen = 295.0 

# high kon
lutfile = joinpath(BASEDIR,"Surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld")

# low kon
# lutfile  = joinpath(BASEDIR,"Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld")

# output control
save_curves = true
visualise   = true

# optimizer parameter ranges (log space except reach)
#logkon_optrange  = (-5.0, -1.25)  # or -2.5 in old file
logkon_optrange = (-3.0, -1.25)
logkoff_optrange = (-4.0, -1.0)
logkonb_optrange = (-3.0, 1.0)
reach_optrange   = (2.0, 35.0)
logCP_optrange   = (1.0, 5.0)

############################ END INPUT ###############################

# the optimizer's parameter ranges (all log space except reach)
optpar_ranges = [logkon_optrange,logkoff_optrange,logkonb_optrange,reach_optrange,logCP_optrange]

# load the surrogate
surrogate = Surrogate(lutfile)

# this should create the directory if it doesn't exist
mkpath(OUTDIR)

# loop through files and do the fitting
allfiles = readdir(RAWDIR)

for file in allfiles
    if occursin(r"^Data_", file) == true
        filename = replace(replace(file,r"^Data_" => "" ), r".csv$" => "")

        # get the aligned SPR data
        aligneddat = get_aligned_data(joinpath(RAWDIR, file), antigenconcen)

        # find the best fit parameters
        bbopt_output = fit_spr_data(surrogate, aligneddat, optpar_ranges)

        if visualise
            figfile = joinpath(OUTDIR, filename * "_fit_curves.png")
            visualisefit(bbopt_output, aligneddat, surrogate.simpars, figfile)
        end
        if save_curves
            curvefile = joinpath(OUTDIR, filename)
            savefit(bbopt_output, aligneddat, surrogate.simpars, curvefile)
        end
    end
end


