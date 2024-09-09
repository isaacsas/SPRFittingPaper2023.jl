# SPRFittingPaper2023.jl API
```@meta
CurrentModule = SPRFittingPaper2023
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
means
means!
vars
sems
SimNumberTerminator
VarianceTerminator
```

## Alignment of SPR Data

```@docs
AlignedData
get_aligned_data
```

## Surrogate

```@docs
SurrogateParams
Surrogate
build_surrogate_serial
save_surrogate
```

## Fitting

```@docs
fit_spr_data
optpars_to_physpars
savefit
visualisefit
```

# Private API Functions

```@docs
surrogate_sprdata_error
update_pars_and_run_spr_sim!
```