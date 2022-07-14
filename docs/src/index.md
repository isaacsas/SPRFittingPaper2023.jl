# SPRFitting.jl

## Running in parallel
1. Save the metadata for the surrogate via modifying and running `examples/make_surrogate_metadata.jl`.
2. Modify and run `examples/make_slices.sh` to batch the surrogate data construction for your cluster.
3. Modify and run `examples/merge_surrogate_slices.jl` to merge the various data files back into the final surrogate.