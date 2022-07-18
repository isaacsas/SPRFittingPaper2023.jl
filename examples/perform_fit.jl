using SPRFitting
using BlackBoxOptim: best_fitness
using LinearAlgebra
using Plots
using OptimizationNLopt   # for fitting monovalent only

############ INPUT ############

# file directories
#BASEDIR = joinpath(@__DIR__, "figures and data")
BASEDIR = "/Users/isaacsas/data/2022-06-07 - FD11A_Data"
RAWDIR = joinpath(BASEDIR, "Aligned")
OUTDIR = joinpath(BASEDIR, "mergetest_mysurrogate_withmono")

# high kon
#lutfile = "/Users/isaacsas/data/surrogates/CombinedLUT_HigherKon_FourParameter_T600_TS150_NG30_Feb4.jld"
lutfile = "/Users/isaacsas/data/surrogates/surrogate_slice_merged.jld"

# low kon
# lutfile  = joinpath(BASEDIR,"Surrogates/CombinedLUT_LowerKon_FourParameter_T600_TS150_NG30_Jan27.jld")

# output control
nfits       = 30   # how many fits to run and then take the minimum over
save_curves = true
visualise   = true
nsims       = 100  # number of simulations to use when plotting

# monovalent fitting, set monovalent_optimizer=nothing if not desired
monovalent_optimizer = NLopt.LD_LBFGS()
lb = [-8.0, -8.0, -8.0]   # upper bounds on parameters in log space (kon,koff,CP)
ub = [8.0, 8.0, 8.0]      # lower bounds on parameters in log space (kon,koff,CP)
ad = Optimization.AutoForwardDiff()   # set to nothing for a derivative-free method
abstol = 1e-8
reltol = 1e-8

# optimizer parameter ranges (log space except reach) for use with surrogate
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
        bbopt_output,best_pars = fit_spr_data(surrogate, aligneddat, optpar_ranges)
        for i in 2:nfits
            println("\nRunning fit: ", i, "/", nfits)
            bbopt_output_new, best_pars_new = fit_spr_data(surrogate, aligneddat, optpar_ranges)
            if best_fitness(bbopt_output_new) < best_fitness(bbopt_output)
                bbopt_output = bbopt_output_new
                best_pars = best_pars_new
            end
        end

        println("Best fit is: ")
        @show best_pars

        # for use with outputting so we don't modify the surrogate's parameters
        simpars = deepcopy(surrogate.simpars)
        simpars.nsims = nsims

        # if we want to include a monovalent fit
        if monovalent_optimizer !== nothing
            # use the bivalent fits as our guess for the monovalent fit
            # [kon, koff, CP]
            u₀ = log10.( [best_pars[1], best_pars[2], best_pars[end]] )
            monofit = monovalent_fit_spr_data(monovalent_optimizer, aligneddat,
                                              simpars.tstop_AtoB, u₀; ad, lb, ub, abstol,
                                              reltol)
        else
            monofit = nothing
        end

        if visualise
            print("saving plot...")
            figfile = joinpath(OUTDIR, filename * "_fit_curves.png")
            visualisefit(bbopt_output, aligneddat, surrogate, simpars, figfile)
            if monofit !== nothing
                figfile = joinpath(OUTDIR, filename * "_fit_curves_monovalent.png")
                visualisefit(monofit, aligneddat, simpars.tstop_AtoB, figfile)
            end
            println("done")
        end
        if save_curves
            print("saving spreadsheet...")
            curvefile = joinpath(OUTDIR, filename)
            savefit(bbopt_output, aligneddat, surrogate, simpars, curvefile; monofit)
            println("done")
        end
    end
end
