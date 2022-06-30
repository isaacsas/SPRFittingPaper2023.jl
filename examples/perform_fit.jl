using SPRFitting
using BlackBoxOptim: best_fitness
using LinearAlgebra
using Plots

############ INPUT ############

# file directories
#BASEDIR = joinpath(@__DIR__, "figures and data")
BASEDIR = "/Users/isaacsas/data/2022-06-07 - FD11A_Data"
RAWDIR = joinpath(BASEDIR, "Aligned")
OUTDIR = joinpath(BASEDIR, "Fitted_Sam_test2")

# high kon
lutfile = "/Users/isaacsas/data/surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld"
#lutfile = joinpath(BASEDIR,"Surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld")

# low kon
# lutfile  = joinpath(BASEDIR,"Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld")

# output control
nfits       = 100   # how many fits to run and then take the minimum over
save_curves = true
visualise   = true
nsims       = 100  # number of simulations to use when plotting

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

for (n,file) in enumerate(allfiles)
    println("\n#####################\n", "Fitting file: ", n, "/", length(allfiles), "\n", "File: ", file, "\n")

    if occursin(r"^Data_", file) == true
        filename = replace(replace(file,r"^Data_" => "" ), r".csv$" => "")

        # get the aligned SPR data
        aligneddat = get_aligned_data(joinpath(RAWDIR, file))

        # find the best fit parameters
        bbopt_output = fit_spr_data(surrogate, aligneddat, optpar_ranges)
        for i in 2:nfits
            println("\nRunning fit: ", i, "/", nfits)
            bbopt_output_new = fit_spr_data(surrogate, aligneddat, optpar_ranges)
            if best_fitness(bbopt_output_new) < best_fitness(bbopt_output)
                bbopt_output = bbopt_output_new
            end
        end

        # for use with outputting so we don't modify the surrogate's parameters
        simpars = deepcopy(surrogate.simpars)
        simpars.nsims = nsims

        if visualise
            print("saving plot...")
            figfile = joinpath(OUTDIR, filename * "_fit_curves.png")
            visualisefit(bbopt_output, aligneddat, surrogate, simpars, figfile)
            println("done")
        end
        if save_curves
            print("saving spreadsheet...")
            curvefile = joinpath(OUTDIR, filename)
            savefit(bbopt_output, aligneddat, surrogate, simpars, curvefile)
            println("done")
        end
    end
end
