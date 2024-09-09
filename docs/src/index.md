# SPRFittingPaper2023.jl

SPRFittingPaper2023.jl provides a library of functions that were used in fitting the particle-based jump process model of [1] to SPR data. The library provides three main components:

1. A forward solver for the particle-based jump process reaction model for
   bivalent antibody-antigen SPR interactions (in both two and three dimesions).
2. Functionality for building the surrogate model approximation the
   particle-based jump process model over a portion of parameter space.
3. Functionality for fitting the surrogate model to SPR data sets to produce
   estimates for the biophysical parameters.

For each of these components we provide a tutorial on their use as part of this documentation. Readers interested in our general methodology should consult [1].

## Installation
To install into your main Julia environment
```julia
using Pkg
Pkg.add(url="https://github.com/isaacsas/SPRFittingPaper2023.jl.git")
```
It is often better to first create a new, clean environment in the directory where
you'll have your fitting script. Start Julia in that directory and then
```julia
using Pkg
Pkg.activate("Environment_Name")
Pkg.add(url="https://github.com/isaacsas/SPRFittingPaper2023.jl.git")
```
You can use the package manager or "Pkg.add" to add any other needed packages to
that environment.


## Bibliography
1. A. Huhn, D. Nissley, ..., C. M. Deane, S. A. Isaacson, and O. Dushek,
   *Analysis of emergent bivalent antibody binding identifies the molecular
   reach as a critical determinant of SARS-CoV-2 neutralisation potency*, in
   review, [available on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.06.556503v2) (2024).