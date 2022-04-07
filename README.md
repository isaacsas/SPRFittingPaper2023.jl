# SPRFitting


## Installation
To install into your main Julia environment
```julia
using Pkg
Pkg.add(url="https://github.com/isaacsas/SPRFitting.jl.git")
```
It is often better to first create a new, clean environment in the directory where
you'll have your fitting script. Start Julia in that directory and then
```julia
using Pkg
Pkg.activate("EnvName")
Pkg.add(url="https://github.com/isaacsas/SPRFitting.jl.git")
```
You can use the package manager or "Pkg.add" to add any other needed packages to
that environment.

## Using the fitting script
Suppose we are in a folder where we've placed and appropriately updated the `examples/perform_fit.jl` script, and run the commands above. We need to also install plots
```julia
Pkg.add("Plots")
```
and then we can run the script as follows:
```julia
include("perform_fit.jl")
```
A similar approach should work for the other scripts in `examples`, modulo that you'll need to `Pkg.add` any libraries that they include via `using` at the top of the script.
