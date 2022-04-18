# SPRFitting.jl API
```@meta
CurrentModule = SPRFitting
```

## Biophysical and Simulation Parameters

```@docs
BioPhysParams
biopars_from_fitting_vec
SimParams
```

## Forward Simulations and Callbacks

```@docs
run_spr_sim!
TotalBoundOutputter
TotalAOutputter
SimNumberTerminator
```

## Alignment of SPR Data

```@docs
AlignedData
get_aligned_data
```

## Surrogate Parameters

```@docs
SurrogateRanges
Surrogate
```

## Fitting 

```@docs
fit_spr_data
bboptpars_to_physpars
savefit
visualisefit
```

