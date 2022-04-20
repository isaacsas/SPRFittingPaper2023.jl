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
VarianceTerminator
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

# Private API Functions

```@docs
surrogate_sprdata_error 
update_pars_and_run_spr_sim!
```