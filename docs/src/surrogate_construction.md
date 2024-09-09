# Surrogate Model Construction
In this tutorial we will illustrate how to construct a (small) surrogate in
serial on any computer, and our workflow for building large surrogates on
clusters (specific to the Grid Engine queuing system based cluster we use).

## Setup
We begin by installing the packages we need in a clean environment:
```julia
using Pkg

# create a new environment for the surrogate constructions
Pkg.activate("surrogate_env") 
Pkg.add(url="https://github.com/isaacsas/SPRFittingPaper2023.jl.git")
Pkg.add("JLD")
```

## Serial Surrogate Construction
We first demonstrate how to build a (small) surrogate in serial (i.e. on a
single CPU core). We start by loading the needed packages:
```@example serialsur
using SPRFittingPaper2023, JLD

# import a useful but non-exported function:
using SPRFittingPaper2023
```

We next define the parameter ranges we wish to tabulate the surrogate over.
Reaction rates are specified via a range in **`log10`** space, i.e.
`log10(kon)`, `log10(koff)`, and `log10(konb)`. The reach is specified in linear
space:
```@example serialsur
logkon_range  = (-3.0, 2.0)
logkoff_range = (-4.0, -1.0)
logkonb_range = (-3.0, 1.5)
reach_range   = (2.0, 35.0)   # in nm
surpars = SurrogateParams(; logkon_range, logkoff_range, logkonb_range, reach_range)
```
We collect these parameters into a [`SurrogateParams`](@ref) structure. Note
that here we will build the surrogate with a fixed internal antibody
concentration of `1.0` nM and antigen concentration of
```@example serialsur
SPRFittingPaper2023.DEFAULT_SIM_ANTIGENCONCEN
```
in units of μM. This is feasible to do because the antibody concentration only
arises in product with $k_{\text{on}}$, so when fitting we can (internally)
interpret the surrogate `logkon` values as representing $\log_{10}(k_{\text{on}}
[\text{Ab}])$. Similarly, as explained in the methods section of [1], using a
fixed antigen concentration is not problematic as we can analytically transform
any fit reach value using the internal antigen concentration to a physical reach
value corresponding to the true experimental antigen concentration. 

Next we specify temporal information for the surrogate:
```@example serialsur
tstop      = 600.0               # simulation end time
tsavelen   = 601                 # number of time points to save (must be integer currently)
tstop_AtoB = 150.0               # time to remove free antibodies
tsave = collect(range(0.0, tstop, tsavelen))  # times to save at
simpars = SimParams(; tstop, tstop_AtoB, tsave)
```
We collect these parameters into a [`SimParams`](@ref) object. Note that it
contains many other parameters for which we typically just use the default value
when building a surrogate (for example, by default distributing antigen
particles uniformly within a cube). The default number of antigen is
```@example serialsur
simpars.N
```

Finally, we specify how many points to tabulate over for each of the four
parameters. Parameters are spaced uniformly in `log10` space for the rates and
linear space for the reach:
```@example serialsur
# here the order is number of points for 
# [logkon, logkoff, logkonb, reach, time]
surrogate_size = (3, 3, 3, 3, tsavelen)
```
Note, this is a very small surrogate, which we would not use in any practical
fitting assay. More typical values are given in our manuscript (often 30-50
points per parameter).

We are now ready to build and save the surrogate
```@example serialsur
outfile = tempname()  # just use a temporary file name
save_surrogate(outfile, surrogate_size, surpars, simpars)
```

Note that the surrogate by default saves curves that correspond to the 
```math
\frac{\text{average number of bound antibodies}}{\text{ the number of antigen in the system}}.
```

## Surrogate Format
The surrogate is stored in a Julia [JLD](https://github.com/JuliaIO/JLD.jl) file. We can see the raw data in the surrogate via the JLD `load`
command:
```@example serialsur
surdata = load(outfile)
```
To access a given field we can say
```@example serialsur
surdata["tstop"]
```
In particular, `surdata["FirstMoment"]` will correspond to the table of solution
curves.

To load the surrogate for use in fitting we instead use
```@example serialsur
sur = Surrogate(outfile)
```
The use of the loaded surrogate for fitting will be illustrated in the next tutorial.


## Parallel Surrogate Construction
Below we explain our workflow for constructing the surrogate via parallel
simulations using the Grid Engine queuing system. The basic workflow and scripts
we link were designed for this system, but should be adaptable to other queue-based 
clusters.

The basic approach is 
1. Save a metadata file with all information needed to build the surrogate.
2. Construct pieces of the surrogate as independent single-core jobs on the cluster.
3. Merge the pieces of the surrogate back together into a single complete surrogate.

The surrogate used in [1] can be downloaded from
[here](https://doi.org/10.6084/m9.figshare.26936854). Below we give the basic
scripts and commands used in its construction (note constructing such a
surrogate generally requires in total 500-2000 hours of cpu time on the Boston
University cluster).

Our workflow is
1. Edit
   ["make\_surrogate\_metadata.jl"](./parallel_surrogate_files/make_surrogate_metadata.jl)
   for your system and then run it in Julia (via
   `include("make_surrogate_metadata.jl")`).
2. Edit the bash scripts ["queue_job.sh"](./parallel_surrogate_files/queue_job.sh)
   and ["run_sim.sh"](./parallel_surrogate_files/run_sim.sh) for your system as
   appropriate and run queue_job.sh to submit the parallel surrogate jobs.
3. After all jobs finish, edit
   ["merge\_surrogate\_slices.jl"](./parallel_surrogate_files/merge_surrogate_slices.jl)
   as appropriate and run it in Julia via `include("merge_surrogate_slices.jl")` to merge the output from each cluster job into one complete surrogate. 

The linked scripts should correspond to those used to construct the surrogate in [1].

## Bibliography
1. A. Huhn, D. Nissley, ..., C. M. Deane, S. A. Isaacson, and O. Dushek,
   *Analysis of emergent bivalent antibody binding identifies the molecular
   reach as a critical determinant of SARS-CoV-2 neutralisation potency*, in
   review, [available on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.06.556503v2) (2024).