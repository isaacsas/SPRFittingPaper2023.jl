# Forward Model Simulation
In this tutorial we will show how to run forward model simulations given a set of parameters and plot the amount of bound antibodies.

## Tutorial Setup
Let's first install the packages we need in a clean environment:
```julia
using Pkg

# create a new environment for the forward model simulation
Pkg.activate("fwd_sim_env") 
Pkg.add(url="https://github.com/isaacsas/SPRFittingPaper2023.jl.git")
Pkg.add("Plots")
Pkg.add("CSV")
```

## Running Forward Simulations
We begin by loading the needed packages:
```@example fwdsim
using SPRFittingPaper2023, Plots, CSV

# import a useful but non-exported function:
using SPRFittingPaper2023: means
```

We next specify the needed model parameters. We will vary the antibody
concentration to generate a sequence of curves to model a series of experiments
that increase the antibody concentration. The parameters we use will be those
that we previously fit to an FD11A-RBD binding experiment.

We will begin using a 3D model in which antigen are randomly distributed within
a cube, our standard model for bivalent SPR experiments. We first specify the
biophysical parameters to use in our simulations and collect them in a
[`BioPhysParams`](@ref) structure:
```@example fwdsim
antigenconcen   = 13.8           # assumed in μM
antibodyconcens = [25.0, 100.0]  # assumed in nM
kon = 5.286e-05                  # assumed in 1 / (nM * sec)
koff = 0.040868575               # assumed in 1 / sec
konb =  0.7801815                # assumed in 1 / sec
reach = 31.89884387              # assumed in nm
CP = 128.569235     # coefficient of proportionality to fit the SPR data
biopars = BioPhysParams(; kon, koff, konb, reach, antigenconcen, CP,
                        antibodyconcen = antibodyconcens[1])
```

!!! note
    Throughout the library, all reported $k_{\text{on}}$ values represent the total bimolecular rate for the binding of free antibodies in solution to an antigen, i.e the $A_i \overset{k_{\text{on}} [\text{Ab}]}{\to} B_i$ reaction, and hence are double the physical $k_{\text{on}}$ defined in our manuscript (where we assume $k_{\text{on}}$ is the rate associated with an individual Fab, i.e. $A_i \overset{2 k_{\text{on}} [\text{Ab}]}{\to} B_i$).

Next we specify simulation parameters. Note, as we gave concentrations for the
antigen and specify the number of particles to use, this fully determines the
domain size (which is calculated for us automatically):
```@example fwdsim
tstop = 600.0        # time in seconds
tstop_AtoB = 150.0    # time to remove free antibodies at
dt = 1.0              # times at which we save the data
N = 1000              # number of antigen particles to use
simpars = SimParams(; antigenconcen = biopars.antigenconcen, tstop, dt, N, 
                      DIM = 3, tstop_AtoB)
```
See the [`SimParams`](@ref) documentation for more on what the various arguments here mean.

Finally, for each antibody concentration we will run `simpars.nsims` forward
simulations and plot the ratio of the average number of bound antibodies to the
number of antigen, scaled by the fitted coefficient of proportionality, `CP`.
Since we set no value for it, by default `simpars.nsims = 1000` are averaged:
```@example fwdsim
plt = plot(; xlabel = "times (sec)", 
             ylabel = "SPR Response")
for abc in antibodyconcens
    biopars.antibodyconcen = abc
    tbo = TotalBoundOutputter(length(simpars.tsave))
    run_spr_sim!(tbo, biopars, simpars)
    bindcnt = means(tbo)
    plot!(plt, simpars.tsave, bindcnt; lw = 2, 
          label = "[Ab] = $(biopars.antibodyconcen) (simulation)")
end

plt
```
Here `means` finalizes calculating the average number of bound antibodies across
the 1000 simulations that were run for each antibody concentration. It then
returns the `CP` scaled average number of antibodies bound to the surface at
each time in `tsave` divided by the number of antigen in the system (i.e. `N =
1000` here). See the [`run_spr_sim!`](@ref) docs and the docs for the
terminators [`SPRFittingPaper2023.VarianceTerminator`](@ref) or
[`SPRFittingPaper2023.SimNumberTerminator`](@ref) for more information on
terminating simulations in relation to the number / accuracy of the desired
averages.

Finally, let's load and plot the corresponding aligned SPR data that we
originally estimated these parameters from to confirm the good fit. 
```@example fwdsim
datadir = joinpath(@__DIR__, "..", "..", "data")
fname = joinpath(datadir, "Data_FC4_10-05-22_Protein07_FD-11A_RBD-13.8_aligned.csv")
ad = get_aligned_data(fname)
for (i, times) in enumerate(ad.times)
    plot!(plt, times, ad.refdata[i]; linestyle = :dash, lw = 2,
          label = "[Ab] = $(ad.antibodyconcens[i]), (data)")
end
plt
```