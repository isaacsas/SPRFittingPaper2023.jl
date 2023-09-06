# SPRFitting - Tools for fitting bivalent antibody SPR assays


## Installation
To install into your main Julia environment
```julia
using Pkg
Pkg.add(url="https://github.com/isaacsas/SPRFittingPaper2023.git")
```
It is often better to first create a new, clean environment in the directory where
you'll have your fitting script. Start Julia in that directory and then
```julia
using Pkg
Pkg.activate("EnvName")
Pkg.add(url="https://github.com/isaacsas/SPRFittingPaper2023.git")
```
You can use the package manager or "Pkg.add" to add any other needed packages to
that environment.

## Creating a Surrogate
The file "examples/make_surrogate.jl" is an example of how to create a surrogate
on one cpu.

The files in "examples/cluster_scripts/" are the typical work flow to create a
surrogate in parallel on a cluster using Grid Engine (qsub). These correspond to
the surrogate used in the manuscript. The workflow is

1. Edit and run "make_surrogate_metadata.jl" to create a file with the shared
   metadata for each parameter set.
2. Edit "queue_job.sh" and "run_sim.sh" to have the appropriate parameters for
   your run. Note "run_sim.h" calls the "examples/make_surrogate_parallel.jl"
   script, so needs to know its location.
3. Edit and run "merge_surrogate_slices.jl" to merge the resulting output files
   into one global surrogate.

## Fitting
Suppose we are in a folder where we've placed and appropriately updated the
"examples/perform_fit.jl" script (make sure to install the packages it uses
too). We can run the script as follows:
```julia
include("perform_fit.jl")
```
SPR data for FD11a that can be fit with a sufficiently resolved surrogate is
available in "data/input_spr_data/*.csv". Note the comments in
"examples/perform_fit.jl" on how the file name must encode the antigen
concentration used for the SPR series within a given csv file.


## Note
When you start a new Julia session in your script folder you'll need to
reactivate the environment you previously created, i.e.
```julia
using Pkg
Pkg.activate("EnvName")
```
where `EnvName` is what you called the environment previously. Running these
commands should cause the environment to be activated and you to be able to run
your scripts again.
