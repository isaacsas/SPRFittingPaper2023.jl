struct Optimiser{Q,R,S,T,U,V}
    """The Optimization.jl method to use."""
    method::Q
    """The AD method to use, nothing by default."""
    ad::R
    """Lower bounds on the parameters; needed by many methods."""
    lb::S
    """Upper bounds on the parameters; needed by many methods."""
    ub::T
    """Named tuple of additional keyword arguments to `OptimizationProblem``."""
    probkwargs::U
    """Name tuple of additional keyword arguments to `solve`. Typically method-specific settings."""
    solverkwargs::V
end

function Optimiser(method; lb = nothing, ub = nothing, ad = nothing,
                           probkwargs = NamedTuple(), solverkwargs = NamedTuple())
    Optimiser(method, ad, lb, ub, probkwargs, solverkwargs)
end


function default_mono_optimiser(lb, ub; kwargs...)
    Optimiser(NLopt.LD_LBFGS(); lb, ub, ad = Optimization.AutoForwardDiff(), kwargs...)
end

function default_bivalent_optimiser(lb, ub; maxiters = 5000, TraceMode = :compact,
                                    TraceInterval = 10.0, probkwargs = NamedTuple(),
                                    solverkwargs = NamedTuple())
    xneskwargs = (; maxiters, TraceMode, TraceInterval)
    solverkwargs = if solverkwargs === nothing
        xneskwargs
    else
        merge(solverkwargs, xneskwargs)
    end
    Optimiser(BBO_xnes(); lb, ub, probkwargs, solverkwargs)
end

@inline scaletoLUT(par, parmin, sz, width) = (par-parmin)*(sz-1)/width + 1

"""
    surrogate_sprdata_error(optpars, surrogate::Surrogate, aligned_data::AlignedData)

This function takes a set of optimization parameters and interpolates a
simulated kinetics curve from the surrogate, returning the ``L^2`` error against
the provided data.

"""
function surrogate_sprdata_error(optpars, pars)
    surrogate = pars.surrogate
    surpars = surrogate.surpars
    sursize = surrogate.surrogate_size

    @unpack times, refdata, antibodyconcens = pars.aligneddat

    # width of each range of parameters in the surrogate
    dq1 = surpars.logkon_range[2] - surpars.logkon_range[1]
    dq2 = surpars.logkoff_range[2] - surpars.logkoff_range[1]
    dq3 = surpars.logkonb_range[2] - surpars.logkonb_range[1]
    dq4 = surpars.reach_range[2] - surpars.reach_range[1]

    # the interpolant assumes a function y = f(x), where we provided y, and x is
    # just an integer for each data point we must therefore rescale from the
    # actual x-value in log parameter space to the integer space used in the
    # interpolant:
    q2 = scaletoLUT(optpars[2], surpars.logkoff_range[1], sursize[2], dq2)
    q3 = scaletoLUT(optpars[3], surpars.logkonb_range[1], sursize[3], dq3)
    q4 = scaletoLUT(optpars[4], surpars.reach_range[1], sursize[4], dq4)

    err = 0.0
    refabc = antibodyconcens[1]
    surrogate_times = surrogate.simpars.tsave
    tspan = last(surrogate_times) - first(surrogate_times)
    @inbounds for (j,abc) in enumerate(antibodyconcens)

        # rescale the on rate to account for changing concentrations of antibody
        logkon = optpars[1] + log10(abc / refabc)
        q1 = scaletoLUT(logkon, surpars.logkon_range[1], sursize[1], dq1)

        sprdata = refdata[j]
        @inbounds for (i,t) in enumerate(times[j])
            s = scaletoLUT(t, surrogate_times[1], length(surrogate_times), tspan)

            # should below be t + 1.0 to scale [0.0,T] to [1,T+1]?
            newerr = 10.0^(optpars[5]) * surrogate.itp(q1,q2,q3,q4,s) - sprdata[i]
            err += newerr * newerr
        end
    end

    # calculate the ℓ₂ error
    return err
end

checkrange(rsur,ropt) = (rsur[1] <= ropt[1] <= ropt[2] <= rsur[2])

function checkranges(optranges, sr::SurrogateParams)
    checkrange(sr.logkon_range, optranges[1]) || error("Optimizer logkon_range not a subset of surrogate logkon_optrange")
    checkrange(sr.logkoff_range, optranges[2]) || error("Optimizer logkoff_range not subset of surrogate logkoff_optrange")
    checkrange(sr.logkonb_range, optranges[3]) || error("Optimizer logkonb_range not subset of surrogate logkonb_optrange")
    checkrange(sr.reach_range, optranges[4]) || error("Optimizer reach_range not subset of surrogate reach_optrange")
    nothing
end

# ensures that the fitting stays within the surrogate by
# shifting down logkon[2] to always be in the desired range
function rescale_logkon_max(logkon::T, sur_logkon, abcs; warn = true) where {T}
    lkshift = log10(maximum(abcs) / abcs[1])
    res = logkon
    if (logkon[2] + lkshift) > sur_logkon[2]
        @set! res[2] = sur_logkon[2] - lkshift
        (res[1] <= res[2]) || error("After the range given by logkon is smaller than needed to handle the provided changes in antibody concentration during fitting.")
        @warn "Decreasing logkon[2] to $(res[2]) to ensure fitting stays within the surrogate."
    end
    res
end

"""
    fit_spr_data(surrogate::Surrogate, aligneddat::AlignedData, searchrange;
                 optimiser = nothing, u₀ = nothing)

Find best fit parameters of the surrogate to the given data.

Notes:
- `searchrange` should be a BlackBoxOptim compatible vector of `Tuple`s of the
  form:
  ```
    searchrange = [logkon_optrange,logkoff_optrange,logkonb_optrange,reach_optrange,logCP_optrange]
  ```
  or
  ```
    searchrange = [logCP_optrange]
  ```
  In the latter case the other parameter ranges are set equal to the range
  within the surrogate.
- Uses `default_bivalent_optimiser` by default if `optimiser = nothing`. Otherwise pass an
  `Optimiser` object.
- `u₀` is an optional guess for the inital search point, methods may or may not use.
- Returns the best fit optimization object and the best fit (bio) parameters as a tuple.
"""
function fit_spr_data(surrogate::Surrogate, aligneddat::AlignedData, searchrange;
                      optimiser = nothing, u₀ = nothing)

    sp = surrogate.surpars
    if length(searchrange) == 1
        sr = [sp.logkon_range, sp.logkoff_range, sp.logkonb_range, sp.reach_range, searchrange[1]]
    else
        sr = deepcopy(searchrange)
    end
    checkranges(sr, sp)

    # adjust logkon[2] to ensure fitting stays within the surrogate
    sr[1] = rescale_logkon_max(sr[1], sp.logkon_range, aligneddat.antibodyconcens)

    lb = map(first, sr)
    ub = map(last, sr)

    if (optimiser === nothing)
        optimiser = default_bivalent_optimiser(lb, ub)
    else
        @set! optimiser.lb = lb
        @set! optimiser.ub = ub
    end

    optfun = if optimiser.ad === nothing
        OptimizationFunction(surrogate_sprdata_error)
    else
        OptimizationFunction(surrogate_sprdata_error, optimiser.ad)
    end

    pars = (; surrogate, aligneddat)
    optprob = OptimizationProblem(optfun, u₀, pars; lb, ub, optimiser.probkwargs...)
    sol = solve(optprob, optimiser.method; optimiser.solverkwargs...)

    # calculate the bestfit biological parameters
    bestpars = optpars_to_physpars(sol.u, aligneddat, surrogate)

    sol, bestpars
end

"""
    optpars_to_physpars(optpars, antibodyconcen, antigenconcen,
                                surrogate_antigenconcen)

    optpars_to_physpars(optpars, aligned_data::AlignedData, surrogate::Surrogate)

Converts parameters vector from BlackBoxOptim to physical parameters, converting
the reach from simulation to physical values.

Notes:
- The reach is converted from the internal parameter value, εᵢ, to the physical
  value, εₑ, via
  ```math
    \\varepsilon_e = \\varepsilon_i \\left(\\frac{[AGC]_i}{[AGC]_e}\\right)^{\\tfrac{1}{3}},
  ```
  where ``[AGC]_i`` is the internal simulator's antigen concentration and
  ``[AGC]_e`` is the concentration used in experiments.
- Assumes these two concentrations have consistent units.
"""
function optpars_to_physpars(optpars, antibodyconcen, antigenconcen,
                               surrogate_antigenconcen)
    kon   = (10.0 ^ optpars[1]) / antibodyconcen  # make bimolecular
    koff  = (10.0 ^ optpars[2])
    konb  = (10.0 ^ optpars[3])
    reach = optpars[4] * cbrt(surrogate_antigenconcen/antigenconcen)
    CP    = (10.0 ^ optpars[5])
    [kon,koff,konb,reach,CP]
end

function optpars_to_physpars(optpars, ad::AlignedData, sur::Surrogate)
    optpars_to_physpars(optpars, ad.antibodyconcens[1], ad.antigenconcen,
                          sur.surpars.antigenconcen)
end


"""
    update_pars_and_run_spr_sim!(outputter, logpars, simpars::SimParams)

Generate the biophysical parameters and run a forward simulation given
parameters from the optimizer.

Arguments:
- outputter = an OutPutter instance for what simulation data to record
- logpars   = vector of the five optimization parameters:
              [logkon,logkoff,logkonb,reach,logCP]
- simpars   = [`SimParams`](@ref) instance, should be consistent with the
              surrogate
"""
function update_pars_and_run_spr_sim!(outputter, logpars, simpars::SimParams)
    # get antigen concentration used in simulations and convert to μM
    antigenconcen = inv_cubic_nm_to_muM(getantigenconcen(simpars))

    # we already account for the antibody concentration in kon, so set it to 1.0
    biopars = biopars_from_fitting_vec(logpars; antigenconcen, antibodyconcen=1.0)

    # reset the outputter
    outputter()

    # run the simulations
    run_spr_sim!(outputter, biopars, simpars)

    nothing
end

"""
    visualisefit(optsol, aligneddat::AlignedData, simpars::SimParams,
                 surrogate::Surrogate, filename=nothing)

Plots fit between data and simulated curves using fitted parameters across a set
of antibody concentrations.

Notes:
- `filename = nothing` if set will cause the graph to be saved.
"""
function visualisefit(optsol, aligneddat::AlignedData, surrogate::Surrogate,
                      simpars::SimParams, filename=nothing, monofit=nothing)
    @unpack times,refdata,antibodyconcens = aligneddat
    params = copy(optsol.u)

    # plot aligned experimental data
    fig1 = plot(; xlabel="time", ylabel="RU", legend=false)
    for (i,sprcurve) in enumerate(refdata)
        plot!(fig1, times[i], sprcurve, color="black")
    end

    # plot simulation data with fit parameters
    ps = copy(params)
    abcref = antibodyconcens[1]
    for (i,abc) in enumerate(antibodyconcens)
        # set kon
        ps[1] = params[1] + log10(abc/abcref)

        # save at the SPR data times
        outputter = TotalBoundOutputter(length(times[i]))
        simpars.tsave = times[i]
        simpars.tstop = last(times[i])

        update_pars_and_run_spr_sim!(outputter, ps, simpars)
        plot!(fig1, times[i], means(outputter))
    end

    (filename !== nothing) && savefig(fig1, filename)

    fig1
end

function visualisefit(monofit, aligneddat::AlignedData, toff::Float64, filename=nothing)
    @unpack times,refdata,antibodyconcens = aligneddat

    # plot aligned experimental data
    fig1 = plot(; xlabel="time", ylabel="RU", legend=false)
    for (i,sprcurve) in enumerate(refdata)
        plot!(fig1, times[i], sprcurve, color="black")
    end

    # plot simulation data with fit parameters
    u = 10.0 .^ monofit.u
    for (i,abc) in enumerate(antibodyconcens)
        ts = times[i]
        totbnd = zeros(length(ts))
        p = (abc,toff)
        for (n,t) in enumerate(ts)
            totbnd[n] = monovalent_total_bound(u, p, t)
        end

        plot!(fig1, times[i], totbnd)
    end

    (filename !== nothing) && savefig(fig1, filename)

    fig1
end


"""
    savefit(optsol, aligneddat::AlignedData, surrogate::Surrogate, simpars::SimParams,
            outfile)

Saves the data, simulated data with fit parameters, and fit parameters in an
XLSX spreadsheet with the given name.
"""
function savefit(optsol, aligneddat::AlignedData, surrogate::Surrogate,
                 simpars::SimParams, outfile; monofit=nothing)
    @unpack times, refdata, antibodyconcens, antigenconcen = aligneddat

    cols_per_time = (monofit === nothing) ? 3 : 4
    savedata = Vector{Vector{Float64}}(undef, cols_per_time*length(antibodyconcens))

    params = copy(optsol.u)
    abcref = antibodyconcens[1]
    toff   = simpars.tstop_AtoB
    ps     = copy(params)
    for (j,abc) in enumerate(antibodyconcens)
        idx = cols_per_time*(j-1) + 1

        savedata[idx] = times[j]
        savedata[idx+1] = refdata[j]

        # times to save model simulations at are same as the SPR data
        outputter = TotalBoundOutputter(length(times[j]))
        simpars.tsave = times[j]
        simpars.tstop = last(times[j])

        # Simulate the model and save the averaged kinetics data
        ps[1] = params[1] + log10(abc/abcref)
        update_pars_and_run_spr_sim!(outputter, ps, simpars)
        savedata[idx+2] = means(outputter)

        # monovalent fit
        if monofit !== nothing
            u = 10.0 .^ monofit.u
            p = (abc, toff)
            totbnd = zeros(length(times[j]))
            for (n,t) in enumerate(times[j])
                totbnd[n] = monovalent_total_bound(u, p, t)
            end
            savedata[idx+3] = totbnd
        end
    end

    # headers for writing the simulation curves
    EH = ["Experimental Data $(antibodyconcens[j]) nM" for j in 1:length(antibodyconcens)]
    MH = ["Bivalent Fit $(antibodyconcens[j]) nM" for j in 1:length(antibodyconcens)]
    MMH = ["Monovalent Fit $(antibodyconcens[j]) nM" for j in 1:length(antibodyconcens)]
    headers = Vector{String}()
    for i in eachindex(antibodyconcens)
        push!(headers, "times")
        push!(headers, EH[i])
        push!(headers, MH[i])
        (monofit !== nothing) && push!(headers, MMH[i])
    end

    # get parameter fits
    bestpars = copy(optsol.u)
    bestfitness = optsol.minimum

    # Write the fits to a spreadsheet
    fname = outfile * "_fit.xlsx"
    XLSX.openxlsx(fname, mode="w") do xf
        # SPR and fit curves
        sheet = xf[1]
        XLSX.rename!(sheet, "SPR and Fit Curves")
        sheet[1,1:length(headers)] = headers
        for (j,sd) in enumerate(savedata)
            nrows = length(sd)
            sheet[2:(nrows+1),j] = sd
        end

        # make separate sheet for fit parameter values
        XLSX.addsheet!(xf, "Fit Parameters")
        sheet = xf[2]

        # internal parameters
        logparnames = ["Best fit parameters (internal):","logkon","logkoff","logkonb","reach","logCP"]
        rows = 1:length(logparnames)
        sheet[rows,1] = logparnames
        sheet[rows,2] = ["", bestpars...]

        # biophysical parameters
        parnames = ["Best fit parameters (physical):","kon","koff","konb","reach","CP"]
        sheet[rows,4] = parnames

        # convert internal simulator concentration from (nm)⁻³ to μM as a consistency check
        simagc = inv_cubic_nm_to_muM(getantigenconcen(simpars))
        @assert isapprox(simagc, surrogate.surpars.antigenconcen, atol=1e-12)

        pars = optpars_to_physpars(bestpars, aligneddat, surrogate)
        sheet[rows,5] = ["", pars...]

        # fitness
        row = last(rows) + 2
        sheet[row,1] = "Bivalent Fitness"
        sheet[row,2] = bestfitness
        if monofit !== nothing
            row += 1
            sheet[row,1] = "Monovalent Fitness"
            sheet[row,2] = monofit.minimum
            sheet[row+2,1] = "Monofit return code:"
            sheet[row+2,2] = string(monofit.retcode)
        end

        # monovalent fit
        if monofit !== nothing
            parnames = ["Monovalent fit parameters (physical):", "kon", "koff", "CP"]
            rows = 1:length(parnames)
            sheet[rows,7] = parnames
            mono_biophys_ps = (10.0 .^ monofit.u)
            sheet[rows,8] = ["", mono_biophys_ps...]
        end
    end

    nothing
end
