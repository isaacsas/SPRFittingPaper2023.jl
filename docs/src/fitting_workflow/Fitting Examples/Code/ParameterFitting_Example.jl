# LIBRARIES
using SPRFitting                    # Custom library for SPR fitting procedures
using Plots
using DataFrames, XLSX, CSV, DelimitedFiles


########################################################################################################################
########################################################################################################################
######################## INPUT ######################## 

experiment_name =  "2024_Example_SPR_curves" #Folder name containing all the SPR files
# all files within the folder will be fit. 
println(experiment_name)

# the algorithm obtains antigen concentration from file name, and antibody concentration from header in SPR file. 
# it is therefore crutial to stick to the filenaming standard. 
# filename: Data_FC[1-4]_[Date]_Protein[No]_[AB name]_Ligand-[AG name]-[AG conc]_aligned.csv
# Example filename: Data_FC2_260224_Protein01_FD-11A_Ligand-RBD-4.44_aligned.csv

# Surrogate
LUTname = # Surrogate file name (LUT = Look up table)
LUTfilename = "surrogate_high_and_low"


# the algorithm obtains lower and upper bounds for bivalent model. 
logCP_optrange   = (1.0, 5.0) # CP range is not specified in LUT and needs to be set extra
optpar_ranges = [logCP_optrange]
# Parameter bounds can also be set manually, if you want to restrict the parameter search space
# optpar_ranges = [(0,1), (1,2), (2,3), (3,4), (4,5)] #[(p1_min,p1_max),(p2_min,p2_max),(p3_min,p3_max),(Lmin,Lmax),(CP_min,CP_max)]


# Fitting 
nfits       = 100   # how many fits to run, the fit with the lowest fitness score is selected.
nsims       = 100  # number of simulations to use when plotting
save_curves = true
visualise   = true


# monovalent fitting, set mono_optimiser=nothing if not desired
# a monovalent model is fit once to determine the quality of the bivalent model fit
lb = [-8.0, -8.0, -8.0]   # lower bounds on parameters in log space (kon,koff,CP)
ub = [8.0, 8.0, 8.0]      # upper bounds on parameters in log space (kon,koff,CP)
mono_optimiser = default_mono_optimiser(lb, ub; solverkwargs = (abstol = 1e-8, reltol = 1e-8))

########################################################################################################################
########################################################################################################################

############ FUNCTIONS ############


function ensure_dir(dir)
   # creats the path if it doesn't exist 
    if isdir(dir) == false
        mkpath(dir)
    end
end

function run_fits(nfits, surrogate, aligneddat, optpar_ranges)
    # function to run fit and find the best fit parameters from nfits runs
    # includes a loop, fit at least twice

    # Arrays to store fitness scores and parameters
    fitness = zeros(nfits)
    physparams = zeros(nfits, 5)
    
    # Perform the first fit
    bbopt_output, best_pars = fit_spr_data(surrogate, aligneddat, optpar_ranges)
    fitness[1] = bbopt_output.minimum
    physparams[1,:] = best_pars
    
    # Iterate over remaining fits
    for i in 2:nfits
        bbopt_output_new, best_pars_new = fit_spr_data(surrogate, aligneddat, optpar_ranges)
        fitness[i] = bbopt_output_new.minimum
        physparams[i,:] = best_pars_new
        if bbopt_output_new.minimum < bbopt_output.minimum
            bbopt_output = bbopt_output_new
            best_pars = best_pars_new
        end
    end
    return fitness, physparams, bbopt_output, best_pars
end


function saveParams_all_fits(params_array, fitness, filename)
    # this saves the returned parameters from all fits, in case you want to see the variability in returned parameters
    SavePara = zeros(Float64,length(fitness),6)
    SavePara[:,1] .= fitness
    SavePara[:,2:end] = params_array

    # Create a DataFrame with the columns
    df = DataFrame(SavePara,:auto)

    Headers = [Symbol("Fitness"), 
               Symbol("kon"), 
               Symbol("koff"),
               Symbol("konb"),
               Symbol("Reach"),
               Symbol("CP")]
    # Rename the headers of each columns
    rename!(df,Headers)

    # Write the DataFrame to an xlsx file
    XLSX.writetable(joinpath(OUTDIR, filename*"_Params.xlsx"),collect(eachcol(df)),names(df), overwrite=true)
end    

function savebestParams(results_paramlist, filenames, expname)
    # this function creates a datafile containing the best paramters for all SPR files in the experiment folder
    df = DataFrame(fnames = filenames, 
               fitness = results_paramlist[:,1],
               kon = results_paramlist[:,2],
               koff = results_paramlist[:,3],
               konb = results_paramlist[:,4],
               reach = results_paramlist[:,5],
               Cp = results_paramlist[:,6]
               )
    CSV.write(joinpath(OUTDIR, expname*"_BestParams.csv"), df)
end

############ RUN SCRIPT ############

# Directories 
BASEDIR = splitdir(dirname(@__FILE__))[1]
EXPDIR = joinpath(BASEDIR,"Experiments", experiment_name)
RAWDIR = joinpath(EXPDIR, "Aligned") 
OUTDIR = joinpath(EXPDIR, "Fitted")

#creating output directory
mkpath(OUTDIR)

#load surrogate
LUT_file = joinpath(BASEDIR,"Surrogates/" * LUTfilename * ".jld")
surrogate = Surrogate(LUT_file)
sps = surrogate.surpars

# Loop through files and do the fitting 
not_hidden(fname::String) = fname[1] != '.' #excludes hidden files starting with .
SPRfiles = filter(not_hidden, readdir(RAWDIR))
println(SPRfiles)

resultcollection = zeros(length(SPRfiles), 6) # array collecting best paramters from all SPR files

for (fx, file) in enumerate(SPRfiles)
    println("\n#####################\n", "Fitting file: ", fx, "/", length(SPRfiles), "\n", "File: ", file, "\n")

    if occursin(r"^Data_", file) == true
        filename = replace(replace(file,r"^Data_" => "" ), r"_aligned.csv$" => "")
        println("Filename", filename, "\n")

        fname = joinpath(RAWDIR, file)
        aligneddat = get_aligned_data(joinpath(RAWDIR, file))

        println("\nRunning fit: ", 1, "/", nfits)

        fitness, physparams, bbopt_output, best_pars = run_fits(nfits, surrogate, aligneddat, optpar_ranges)

        println("Best fit is: ")
        @show best_pars

        # save best result for all files
        resultcollection[fx,1] = bbopt_output.minimum
        resultcollection[fx,2:end] = best_pars
 
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
            visualisefit(bbopt_output, aligneddat, surrogate, simpars, figfile)
            if monofit !== nothing
                figfile = joinpath(OUTDIR, filename * "_fit_curves_monovalent.png")
                visualisefit(monofit, aligneddat, simpars.tstop_AtoB, figfile)
            end
            println("done")
        end
        if save_curves
            print("saving spreadsheet and parameters...")
            curvefile = joinpath(OUTDIR, filename)
            savefit(bbopt_output, aligneddat, surrogate, simpars, curvefile; monofit)
            saveParams_all_fits(physparams, fitness, filename)
            println("done")
        end

         
    end
end

#saving the best parameters of all SPR files in a final data table
savebestParams(resultcollection, SPRfiles, experiment_name)

