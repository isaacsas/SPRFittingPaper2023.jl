var documenterSearchIndex = {"docs":
[{"location":"forward_simulation/#Forward-Model-Simulation","page":"Forward Model","title":"Forward Model Simulation","text":"","category":"section"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"In this tutorial we will show how to run forward model simulations given a set of parameters and plot the amount of bound antibodies.","category":"page"},{"location":"forward_simulation/#Tutorial-Setup","page":"Forward Model","title":"Tutorial Setup","text":"","category":"section"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"Let's first install the packages we need in a clean environment:","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"using Pkg\n\n# create a new environment for the forward model simulation\nPkg.activate(\"fwd_sim_env\") \nPkg.add(url=\"https://github.com/isaacsas/SPRFittingPaper2023.jl.git\")\nPkg.add(\"Plots\")\nPkg.add(\"CSV\"\"","category":"page"},{"location":"forward_simulation/#Running-Forward-Simulations","page":"Forward Model","title":"Running Forward Simulations","text":"","category":"section"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"We begin by loading the needed packages:","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"using SPRFittingPaper2023, Plots, CSV\n\n# import a useful but non-exported function:\nusing SPRFittingPaper2023: means","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"We next specify the needed model parameters. We will vary the antibody concentration to generate a sequence of curves to model a series of experiments that increase the antibody concentration. The parameters we use will be those that we previously fit to an FD11A-RBD binding experiment.","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"We will begin using a 3D model in which antigen are randomly distributed within a cube, our standard model for bivalent SPR experiments. We first specify the biophysical parameters to use in our simulations and collect them in a BioPhysParams structure:","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"antigenconcen   = 13.8           # assumed in μM\nantibodyconcens = [25.0, 100.0]  # assumed in nM\nkon = 5.286e-05                  # assumed in 1 / (nM * sec)\nkoff = 0.040868575               # assumed in 1 / sec\nkonb =  0.7801815                # assumed in 1 / sec\nreach = 31.89884387              # assumed in nm\nCP = 128.569235     # coefficient of proportionality to fit the SPR data\nbiopars = BioPhysParams(; kon, koff, konb, reach, antigenconcen, CP,\n                        antibodyconcen = antibodyconcens[1])","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"note: Note\nThroughout the library, all reported k_texton values represent the total bimolecular rate for the binding of free antibodies in solution to an antigen, i.e the A_i oversetk_texton textAbto B_i reaction, and hence are double the physical k_texton defined in our manuscript (where we assume k_texton is the rate associated with an individual Fab, i.e. A_i overset2 k_texton textAbto B_i).","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"Next we specify simulation parameters. Note, as we gave concentrations for the antigen and specify the number of particles to use, this fully determines the domain size (which is calculated for us automatically):","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"tstop = 600.0        # time in seconds\ntstop_AtoB = 150.0    # time to remove free antibodies at\ndt = 1.0              # times at which we save the data\nN = 1000              # number of antigen particles to use\nsimpars = SimParams(; antigenconcen = biopars.antigenconcen, tstop, dt, N, \n                      DIM = 3, tstop_AtoB)","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"See the SimParams documentation for more on what the various arguments here mean.","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"Finally, for each antibody concentration we will run one forward simulation and plot the ratio of the number of bound antibodies to the number of antigen, scaled by the fitted coefficient of proportionality, CP:","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"plt = plot(; xlabel = \"times (sec)\", \n             ylabel = \"SPR Response\")\nfor abc in antibodyconcens\n    biopars.antibodyconcen = abc\n    tbo = TotalBoundOutputter(length(simpars.tsave))\n    run_spr_sim!(tbo, biopars, simpars)\n    bindcnt = means(tbo)\n    plot!(plt, simpars.tsave, bindcnt; lw = 2, \n          label = \"[Ab] = $(biopars.antibodyconcen) (simulation)\")\nend\n\nplt","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"Note that here we only ran one simulation for each parameter set, and hence means just returns the CP scaled number of antibodies bound to the surface at each time in tsave divided by the number of antigen in the system (i.e. N = 1000 here). If we had wanted to run and average multiple samples we could have passed run_spr_sim! a terminator such as a SPRFittingPaper2023.VarianceTerminator or SPRFittingPaper2023.SimNumberTerminator.","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"Finally, let's load and plot the corresponding aligned SPR data that we originally estimated these parameters from to confirm the good fit. ","category":"page"},{"location":"forward_simulation/","page":"Forward Model","title":"Forward Model","text":"datadir = joinpath(@__DIR__, \"..\", \"..\", \"data\")\nfname = joinpath(datadir, \"Data_FC4_10-05-22_Protein07_FD-11A_RBD-13.8_aligned.csv\")\nad = get_aligned_data(fname)\nfor (i, times) in enumerate(ad.times)\n    plot!(plt, times, ad.refdata[i]; linestyle = :dash, lw = 2,\n          label = \"[Ab] = $(ad.antibodyconcens[i]), (data)\")\nend\nplt","category":"page"},{"location":"surrogate_construction/#Surrogate-Model-Construction","page":"Surrogate Construction","title":"Surrogate Model Construction","text":"","category":"section"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"In this tutorial we will illustrate how to construct a (small) surrogate in serial on any computer, and our workflow for building large surrogates on clusters (specific to the Grid Engine queuing system based cluster we use).","category":"page"},{"location":"surrogate_construction/#Setup","page":"Surrogate Construction","title":"Setup","text":"","category":"section"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"We begin by installing the packages we need in a clean environment:","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"using Pkg\n\n# create a new environment for the forward model simulation\nPkg.activate(\"fwd_sim_env\") \nPkg.add(url=\"https://github.com/isaacsas/SPRFittingPaper2023.jl.git\")\nPkg.add(\"JLD\")","category":"page"},{"location":"surrogate_construction/#Serial-Surrogate-Construction","page":"Surrogate Construction","title":"Serial Surrogate Construction","text":"","category":"section"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"We first demonstrate how to build a (small) surrogate in serial (i.e. on a single CPU core). We start by loading the needed packages:","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"using SPRFittingPaper2023, JLD\n\n# import a useful but non-exported function:\nusing SPRFittingPaper2023","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"We next define the parameter ranges we wish to tabulate the surrogate over. Reaction rates are specified via a range in log10 space, i.e. log10(kon), log10(koff), and log10(konb). The reach is specified in linear space:","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"logkon_range  = (-3.0, 2.0)\nlogkoff_range = (-4.0, -1.0)\nlogkonb_range = (-3.0, 1.5)\nreach_range   = (2.0, 35.0)   # in nm\nsurpars = SurrogateParams(; logkon_range, logkoff_range, logkonb_range, reach_range)","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"We collect these parameters into a SurrogateParams structure. Note that here we will build the surrogate with a fixed internal antibody concentration of 1.0 nM and antigen concentration of","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"SPRFittingPaper2023.DEFAULT_SIM_ANTIGENCONCEN","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"in units of μM. This is feasible to do because the antibody concentration only arises in product with k_texton, so when fitting we can (internally) interpret the surrogate logkon values as representing log_10(k_texton textAb). Similarly, as explained in the methods section of [1], using a fixed antigen concentration is not problematic as we can analytically transform any fit reach value using the internal antigen concentration to a physical reach value corresponding to the true experimental antigen concentration. ","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Next we specify temporal information for the surrogate:","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"tstop      = 600.0               # simulation end time\ntsavelen   = 601                 # number of time points to save (must be integer currently)\ntstop_AtoB = 150.0               # time to remove free antibodies\ntsave = collect(range(0.0, tstop, tsavelen))  # times to save at\nsimpars = SimParams(; tstop, tstop_AtoB, tsave)","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"We collect these parameters into a SimParams object. Note that it contains many other parameters for which we typically just use the default value when building a surrogate (for example, by default distributing antigen particles uniformly within a cube). The default number of antigen is","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"simpars.N","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Finally, we specify how many points to tabulate over for each of the four parameters. Parameters are spaced uniformly in log10 space for the rates and linear space for the reach:","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"# here the order is number of points for \n# [logkon, logkoff, logkonb, reach, time]\nsurrogate_size = (3, 3, 3, 3, tsavelen)","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Note, this is a very small surrogate, which we would not use in any practical fitting assay. More typical values are given in our manuscript (often 30-50 points per parameter).","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"We are now ready to build and save the surrogate","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"outfile = tempname()  # just use a temporary file name\nsave_surrogate(outfile, surrogate_size, surpars, simpars)","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Note that the surrogate by default saves curves that correspond to the ","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"fractextaverage number of bound antibodiestext the number of antigen in the system","category":"page"},{"location":"surrogate_construction/#Surrogate-Format","page":"Surrogate Construction","title":"Surrogate Format","text":"","category":"section"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"The surrogate is stored in a Julia JLD file. We can see the raw data in the surrogate via the JLD load command:","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"surdata = load(outfile)","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"To access a given field we can say","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"surdata[\"tstop\"]","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"In particular, surdata[\"FirstMoment\"] will correspond to the table of solution curves.","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"To load the surrogate for use in fitting we instead use","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"sur = Surrogate(outfile)","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"The use of the loaded surrogate for fitting will be illustrated in the next tutorial.","category":"page"},{"location":"surrogate_construction/#Parallel-Surrogate-Construction","page":"Surrogate Construction","title":"Parallel Surrogate Construction","text":"","category":"section"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Below we explain our workflow for constructing the surrogate via parallel simulations using the Grid Engine queuing system. The basic workflow and scripts we link were designed for this system, but should be adaptable to other queue-based  clusters.","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"The basic approach is ","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Save a metadata file with all information needed to build the surrogate.\nConstruct pieces of the surrogate as independent single-core jobs on the cluster.\nMerge the pieces of the surrogate back together into a single complete surrogate.","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"The surrogate used in [1] can be downloaded from here. Below we give the basic scripts and commands used in its construction (note constructing such a surrogate generally requires in total 500-2000 hours of cpu time on the Boston University cluster).","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Our workflow is","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"Edit \"make_surrogate_metadata.jl\" for your system and then run it in Julia (via include(\"make_surrogate_metadata.jl\")).\nEdit the bash scripts \"queue_job.sh\" and \"run_sim.sh\" for your system as appropriate and run queue_job.sh to submit the parallel surrogate jobs.\nAfter all jobs finish, edit \"merge_surrogate_slices.jl\" as appropriate and run it in Julia via include(\"merge_surrogate_slices.jl\") to merge the output from each cluster job into one complete surrogate. ","category":"page"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"The linked scripts should correspond to those used to construct the surrogate in [1].","category":"page"},{"location":"surrogate_construction/#Bibliography","page":"Surrogate Construction","title":"Bibliography","text":"","category":"section"},{"location":"surrogate_construction/","page":"Surrogate Construction","title":"Surrogate Construction","text":"A. Huhn, D. Nissley, ..., C. M. Deane, S. A. Isaacson, and O. Dushek, Analysis of emergent bivalent antibody binding identifies the molecular reach as a critical determinant of SARS-CoV-2 neutralisation potency, in review, available on bioRxiv (2024).","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.jl-API","page":"API","title":"SPRFittingPaper2023.jl API","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"CurrentModule = SPRFittingPaper2023","category":"page"},{"location":"SPRFitting_api/#Biophysical-and-Simulation-Parameters","page":"API","title":"Biophysical and Simulation Parameters","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"BioPhysParams\nbiopars_from_fitting_vec\nSimParams","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.BioPhysParams","page":"API","title":"SPRFittingPaper2023.BioPhysParams","text":"mutable struct BioPhysParams{T<:Number}\n\nBiophysical parameters to use in forward simulations.\n\nFields\n\nkon: A –> B rate, units of (μM s)⁻¹\nkoff: B –> A rate, units of s⁻¹\nkonb: A+B –> C rate, units of s⁻¹\nreach: Reach of A+B –> C reaction, units of nm\nCP: CP factor (default = 1.0)\nantigenconcen: Concentration of antigen (default = DEFAULT_SIM_ANTIGENCONCEN μM)\nantibodyconcen: Concentration of antibodies, should be consistent with kon's units (default = 1.0 μM)\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.biopars_from_fitting_vec","page":"API","title":"SPRFittingPaper2023.biopars_from_fitting_vec","text":"biopars_from_fitting_vec(p; antibodyconcen=1.0, antigenconcen=1.0)\n\nGenerate a BioPhysParams struct from a parameter vector from fitting.\n\nNotes: \n\nAssumed that    p = log_10(textkon) log_10(textkoff) log_10(textkon_textb) textreach log_10(textCP)\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.SimParams","page":"API","title":"SPRFittingPaper2023.SimParams","text":"mutable struct SimParams{T<:Number, DIM}\n\nBiophysical parameters to use in forward simulations.\n\nFields\n\nN: Number of particles (default = 1000)\ntstop: Time to end simulations at (default = 600.0)\ntstop_AtoB: Time to turn off A –> B reaction, i.e. time the antibody bath is removed. (default = Inf)\ntsave: Times to save data at (default = nothing)\nL: Domain Length\ninitlocs: Initial Particle Positions (default = uniform in -L2L2^textDIM)\nresample_initlocs: Resample initlocs every simulation (default = true)\nnsims: Number of simulations in runsim (to control sampling error) (default = 1000)\n\nKeyword Arguments:\n\nAll fields have a corresponding kwarg.\nantigenconcen = DEFAULT_SIM_ANTIGENCONCEN, the default antigenconcen in units of (nm)⁻³ (unless using convert_agc_units=false).\nDIM = 3 the dimension of the underlying space (i.e. 2 or 3). \nconvert_agc_units = true, set to false to disable the conversion of the antigen units from assumed units of μM to (nm)⁻³, see below. In this case L is calculated using the antigen concentration value directly. \nresample_initlocs = true, if set to false initlocs will be constant, simply reusing the initial value sampled in the SimParams constructor.\n\nNotes:\n\nUses the antigen concentration to determine L, and so needs to be updated if this concentration changes. The antigen concentration is first converted from assumed units of  μM to units of (nm)⁻³. Then \n  L = left(fracNtextantigen text in (nm)^-3right)^tfrac1DIM\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#Forward-Simulations-and-Callbacks","page":"API","title":"Forward Simulations and Callbacks","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"run_spr_sim!\nTotalBoundOutputter\nTotalAOutputter\nmeans\nmeans!\nvars\nsems\nSimNumberTerminator\nVarianceTerminator","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.run_spr_sim!","page":"API","title":"SPRFittingPaper2023.run_spr_sim!","text":"run_spr_sim!(outputter, biopars::BioPhysParams, simpars::SimParams,\n             terminator = SimNumberTerminator())\n\nRuns a set of SPR simulations via the particle model for a fixed set of physical parameters.\n\nNotes:\n\nUse the passed BioPhysParams and SimParams to set simulation parameters (including number of simulations).\nUse the outputter to specify a callback controlling what data is saved.\nUse the terminator to specify a callback controlling when the simulation is stopped.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.TotalBoundOutputter","page":"API","title":"SPRFittingPaper2023.TotalBoundOutputter","text":"struct TotalBoundOutputter{T}\n\nCallback that saves the amount of bound antibodies (i.e. monovalently + bivalently bound) at each save time.\n\nFields\n\nbindcnt\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.TotalAOutputter","page":"API","title":"SPRFittingPaper2023.TotalAOutputter","text":"struct TotalAOutputter{T}\n\nCallback that saves the amount of unbound antigen at each save time.\n\nFields\n\nbindcnt\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.means","page":"API","title":"SPRFittingPaper2023.means","text":"means(o::TotalBoundOutputter)\n\nVector of means.\n\n\n\n\n\nmeans(o::TotalAOutputter)\n\nVector of means.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.means!","page":"API","title":"SPRFittingPaper2023.means!","text":"means!(m, o::TotalBoundOutputter)\n\nIn-place vector of means.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.vars","page":"API","title":"SPRFittingPaper2023.vars","text":"vars(o::TotalBoundOutputter)\n\nVector of variances.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.sems","page":"API","title":"SPRFittingPaper2023.sems","text":"means(o::TotalBoundOutputter)\n\nVector of standard errors.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.SimNumberTerminator","page":"API","title":"SPRFittingPaper2023.SimNumberTerminator","text":"mutable struct SimNumberTerminator\n\nCallback that stops simulating when the desired number of simulations is reached.\n\nFields\n\nnum_completed_sims: How many simulations have been completed\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.VarianceTerminator","page":"API","title":"SPRFittingPaper2023.VarianceTerminator","text":"mutable struct VarianceTerminator\n\nCallback that stops simulating when the variance in bound antibodies becomes sufficiently small or a maximum number of simulations is reached.\n\nFields\n\nssetol: The tolerance below which to stop simulating (default = .01).\nnotdone: True if should keep iterating.\nminsims: Minimum number of sims to run (default = 15).\nmaxsims: Maximum number of sims to run (default = 250).\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#Alignment-of-SPR-Data","page":"API","title":"Alignment of SPR Data","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"AlignedData\nget_aligned_data","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.AlignedData","page":"API","title":"SPRFittingPaper2023.AlignedData","text":"struct AlignedData\n\nThe aligned data from a set of SPR experiments varying antibody concentrations.\n\nFields\n\ntimes: times[i] gives the SPR data timepoints for antibodyconcen[i]\nrefdata: refdata[i] gives the SPR simulation data for antibodyconcens[i]\nantibodyconcens: Antibody concentrations for the SPR data set in μM\nantigenconcen: Antigen concentration for the SPR data set in μM\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.get_aligned_data","page":"API","title":"SPRFittingPaper2023.get_aligned_data","text":"get_aligned_data(fname, antigenconcen)\n\nRead the given CSV file representing aligned SPR data. Returns an AlignedData structure storing the given SPR experiment data.\n\nNotes:\n\nantigenconcen = nothing the concentration of antigen used in the experiments in units of μM. If not set assumes the CSV filename has the form \"\"text-AGCCONCEN_aligned.csv\"\", and parses the number antigen concentration from AGCCONCEN.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#Surrogate","page":"API","title":"Surrogate","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"SurrogateParams\nSurrogate\nbuild_surrogate_serial\nsave_surrogate","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.SurrogateParams","page":"API","title":"SPRFittingPaper2023.SurrogateParams","text":"struct SurrogateParams{T<:Number}\n\nBiophysical parameters used in the surrogate.\n\nFields\n\nlogkon_range: log₁₀ space range of kon, first order rate in surrogate (units of s⁻¹).\nlogkoff_range: log₁₀ space range of koff (units of s⁻¹).\nlogkonb_range: log₁₀ space range of konb, Doi association rate in surrogate (units of s⁻¹).\nreach_range: linear range of reach values (units of nm)\nantigenconcen: Surrogate's internal antigen concentration in μM (default is DEFAULT_SIM_ANTIGENCONCEN).\nantibodyconcen: Surrogate's internal antibody concentration in μM (default is 1.0).\nCP: Surrogate's internal CP value (default is 1.0)\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.Surrogate","page":"API","title":"SPRFittingPaper2023.Surrogate","text":"mutable struct Surrogate{S, T, U, V, W}\n\nSurrogate parameters and data.\n\nFields\n\nsurrogate_size: Number of points for each coordinate within the surrogate lookup table.\nitp: Interpolation of the surrogate lookup table.\nsurpars: SurrogateParams used in the surrogate.\nsimpars: Simulation parameters used in the surrogate.\n\nArguments (one of the following two):\n\nlutfile = the name of the file storing a surrogate to load.\nlutdata = AbstractArray representing the raw data points to build the surrogate from.\n\nNotes:\n\nAll fields can also be passed as a keyword arg.\n\n\n\n\n\n","category":"type"},{"location":"SPRFitting_api/#SPRFittingPaper2023.build_surrogate_serial","page":"API","title":"SPRFittingPaper2023.build_surrogate_serial","text":"build_surrogate_serial(surrogate_size::Tuple, surpars::SurrogateParams, simpars::SimParams;\n                       terminator=VarianceTerminator())\n\nConstruct the data table array for a new surrogate by varying kon, koff, konb and reach.\n\nArguments:\n\nsurrogate_size = a Tuple with (nkon,nkoff,nkonb,nreach) points to use.\nsurpars = the physical parameters to use in the surrogate. It is strongly recommended to not change the default values of antigenconcen, antibodyconcen, or CP unless you really know what you are doing –  these default values are implicitly assumed in other places.\nsimpars = the simulation parameters to use (number of particles, simulations, domain size, etc).\n\nKeyword Arguments:\n\nterminator, can be used to alter how the number of samples for each parameter set is determined. See VarianceTerminator for the default values.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.save_surrogate","page":"API","title":"SPRFittingPaper2023.save_surrogate","text":"save_surrogate(filename, surrogate_size, surpars::SurrogateParams, simpars::SimParams)\n\nCreate a JLD file with the given surrogate metadata.\n\nArguments:\n\nsurrogate_size = a Tuple with (nkon,nkoff,nkonb,nreach) points to use.\nsurpars = the physical parameters to use in the surrogate. It is strongly recommended to not change the default values of antigenconcen, antibodyconcen, or CP unless you really know what you are doing –  these default values are implicitly assumed in other places.\nsimpars = the simulation parameters to use (number of particles, simulations, domain size, etc).\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#Fitting","page":"API","title":"Fitting","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"fit_spr_data\noptpars_to_physpars\nsavefit\nvisualisefit","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.fit_spr_data","page":"API","title":"SPRFittingPaper2023.fit_spr_data","text":"fit_spr_data(surrogate::Surrogate, aligneddat::AlignedData, searchrange;\n             optimiser = nothing, u₀ = nothing)\n\nFind best fit parameters of the surrogate to the given data.\n\nNotes:\n\nsearchrange should be a BlackBoxOptim compatible vector of Tuples of the form:\n  searchrange = [logkon_optrange,logkoff_optrange,logkonb_optrange,reach_optrange,logCP_optrange]\nor\n  searchrange = [logCP_optrange]\nIn the latter case the other parameter ranges are set equal to the range within the surrogate.\nUses default_bivalent_optimiser by default if optimiser = nothing. Otherwise pass an Optimiser object.\nu₀ is an optional guess for the inital search point, methods may or may not use.\nReturns the best fit optimization object and the best fit (bio) parameters as a tuple.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.optpars_to_physpars","page":"API","title":"SPRFittingPaper2023.optpars_to_physpars","text":"optpars_to_physpars(optpars, antibodyconcen, antigenconcen,\n                            surrogate_antigenconcen)\n\noptpars_to_physpars(optpars, aligned_data::AlignedData, surrogate::Surrogate)\n\nConverts parameters vector from BlackBoxOptim to physical parameters, converting the reach from simulation to physical values.\n\nNotes:\n\nThe reach is converted from the internal parameter value, εᵢ, to the physical value, εₑ, via\n  varepsilon_e = varepsilon_i left(fracAGC_iAGC_eright)^tfrac13\nwhere AGC_i is the internal simulator's antigen concentration and AGC_e is the concentration used in experiments.\nAssumes these two concentrations have consistent units.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.savefit","page":"API","title":"SPRFittingPaper2023.savefit","text":"savefit(optsol, aligneddat::AlignedData, surrogate::Surrogate, simpars::SimParams,\n        outfile)\n\nSaves the data, simulated data with fit parameters, and fit parameters in an XLSX spreadsheet with the given name.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.visualisefit","page":"API","title":"SPRFittingPaper2023.visualisefit","text":"visualisefit(optsol, aligneddat::AlignedData, simpars::SimParams,\n             surrogate::Surrogate, filename=nothing)\n\nPlots fit between data and simulated curves using fitted parameters across a set of antibody concentrations.\n\nNotes:\n\nfilename = nothing if set will cause the graph to be saved.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#Private-API-Functions","page":"API","title":"Private API Functions","text":"","category":"section"},{"location":"SPRFitting_api/","page":"API","title":"API","text":"surrogate_sprdata_error\nupdate_pars_and_run_spr_sim!","category":"page"},{"location":"SPRFitting_api/#SPRFittingPaper2023.surrogate_sprdata_error","page":"API","title":"SPRFittingPaper2023.surrogate_sprdata_error","text":"surrogate_sprdata_error(optpars, surrogate::Surrogate, aligned_data::AlignedData)\n\nThis function takes a set of optimization parameters and interpolates a simulated kinetics curve from the surrogate, returning the L^2 error against the provided data.\n\n\n\n\n\n","category":"function"},{"location":"SPRFitting_api/#SPRFittingPaper2023.update_pars_and_run_spr_sim!","page":"API","title":"SPRFittingPaper2023.update_pars_and_run_spr_sim!","text":"update_pars_and_run_spr_sim!(outputter, logpars, simpars::SimParams)\n\nGenerate the biophysical parameters and run a forward simulation given parameters from the optimizer.\n\nArguments:\n\noutputter = an OutPutter instance for what simulation data to record\nlogpars   = vector of the five optimization parameters:             [logkon,logkoff,logkonb,reach,logCP]\nsimpars   = SimParams instance, should be consistent with the             surrogate\n\n\n\n\n\n","category":"function"},{"location":"#SPRFittingPaper2023.jl","page":"Home","title":"SPRFittingPaper2023.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SPRFittingPaper2023.jl provides a library of functions that were used in fitting the particle-based jump process model of [1] to SPR data. The library provides three main components:","category":"page"},{"location":"","page":"Home","title":"Home","text":"A forward solver for the particle-based jump process reaction model for bivalent antibody-antigen SPR interactions (in both two and three dimesions).\nFunctionality for building the surrogate model approximation the particle-based jump process model over a portion of parameter space.\nFunctionality for fitting the surrogate model to SPR data sets to produce estimates for the biophysical parameters.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For each of these components we provide a tutorial on their use as part of this documentation. Readers interested in our general methodology should consult [1].","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install into your main Julia environment","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(url=\"https://github.com/isaacsas/SPRFittingPaper2023.jl.git\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is often better to first create a new, clean environment in the directory where you'll have your fitting script. Start Julia in that directory and then","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.activate(\"Environment_Name\")\nPkg.add(url=\"https://github.com/isaacsas/SPRFittingPaper2023.jl.git\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can use the package manager or \"Pkg.add\" to add any other needed packages to that environment.","category":"page"},{"location":"#Bibliography","page":"Home","title":"Bibliography","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A. Huhn, D. Nissley, ..., C. M. Deane, S. A. Isaacson, and O. Dushek, Analysis of emergent bivalent antibody binding identifies the molecular reach as a critical determinant of SARS-CoV-2 neutralisation potency, in review, available on bioRxiv (2024).","category":"page"}]
}
