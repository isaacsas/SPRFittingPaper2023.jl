using SPRFittingPaper2023
using BlackBoxOptim: best_fitness
using LinearAlgebra
using Plots

############ INPUT ############

# file directories
BASEDIR = joinpath("..", "data")

# folder of input SPR data to fit
# we fit each input csv in the folder
# assumes fname has the form "anytext-ANTIGENCONCENVALUE_aligned.csv"
RAWDIR = joinpath(BASEDIR, "input_spr_data")

# location of output fitted parameter spreadsheets
OUTDIR = joinpath(BASEDIR, "test_output")

# surrogate file
lutfile = joinpath(BASEDIR, "surrogate_slice_merged.jld")

# output control
nfits       = 100   # how many fits to run and then take the minimum over
save_curves = true
visualise   = true
nsims       = 100   # number of simulations to use when plotting

# monovalent fitting, set mono_optimiser=nothing if not desired
lb = [-8.0, -8.0, -8.0]   # lower bounds on parameters in log space (kon,koff,CP)
ub = [8.0, 8.0, 8.0]      # upper bounds on parameters in log space (kon,koff,CP)
mono_optimiser = default_mono_optimiser(lb, ub; solverkwargs = (abstol = 1e-8, reltol = 1e-8))

# optimizer parameter ranges (log space except reach) for use with surrogate
logCP_optrange   = (1.0, 5.0)

############################ END INPUT ###############################

# load the surrogate
surrogate = Surrogate(lutfile)
sps = surrogate.surpars

# we have to set the logCP range here as it is not in the surrogate
optpar_ranges = [logCP_optrange]

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
        optsol, best_pars = fit_spr_data(surrogate, aligneddat, optpar_ranges)
        for i in 2:nfits
            println("\nRunning fit: ", i, "/", nfits)
            optsol_new, best_pars_new = fit_spr_data(surrogate, aligneddat, optpar_ranges)
            if optsol_new.minimum < optsol.minimum
                optsol = optsol_new
                best_pars = best_pars_new
            end
        end

        println("Best fit is: ")
        @show best_pars

        # for use with outputting so we don't modify the surrogate's parameters
        simpars = deepcopy(surrogate.simpars)
        simpars.nsims = nsims

        # if we want to include a monovalent fit
        if mono_optimiser !== nothing
            # use the bivalent fits as our guess for the monovalent fit
            # [kon, koff, CP]
            u₀ = log10.( [best_pars[1], best_pars[2], best_pars[end]] )
            monofit = monovalent_fit_spr_data(mono_optimiser, aligneddat, simpars.tstop_AtoB, u₀)
        else
            monofit = nothing
        end

        if visualise
            print("saving plot...")
            figfile = joinpath(OUTDIR, filename * "_fit_curves.png")
            visualisefit(optsol, aligneddat, surrogate, simpars, figfile)
            if monofit !== nothing
                figfile = joinpath(OUTDIR, filename * "_fit_curves_monovalent.png")
                visualisefit(monofit, aligneddat, simpars.tstop_AtoB, figfile)
            end
            println("done")
        end
        if save_curves
            print("saving spreadsheet...")
            curvefile = joinpath(OUTDIR, filename)
            savefit(optsol, aligneddat, surrogate, simpars, curvefile; monofit)
            println("done")
        end
    end
end
