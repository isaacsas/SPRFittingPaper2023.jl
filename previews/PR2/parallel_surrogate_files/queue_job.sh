#! /bin/bash

# the path to where you have installed SPRFittingPaper2023
# In Julia you can get this via calling pkgdir(SPRFittingPaper2023)
PKGPATH="/project/fpkmc3d/.julia_sai/dev/SPRFitting"

# The output name for the metadata file you created
METADATAFILE="/project/fpkmc3d/surrogates/rebuild_lut_highandlow/surrogate_metadata.jld"

# the base name to append to each surrogate portion (i.e. what each cluster CPU builds)
OUTFILEBASENAME="/project/fpkmc3d/surrogates/rebuild_lut_highandlow/surrogate_slice"

# the julia script we run to build one subset of the surrogate
JULIASCRIPT="${PKGPATH}/examples/make_surrogate_parallel.jl"

# the total number of parameter sets we iterate over, in this case
# 42*30*30*30 
NUMPARAMSETS=1134000

# the number of independent jobs to create, each job handles NUMPARAMSETS / NCPUs parameter sets
NCPUS=900

let IDXINCREMENT=$NUMPARAMSETS/$NCPUS
let dn=$IDXINCREMENT-1
let idx=1
for (( n=1;n<=$NUMPARAMSETS;n+=$IDXINCREMENT)); do
    let last=$n+$dn
    qsub run_sim.sh "$JULIASCRIPT" "$PKGPATH" "$METADATAFILE" "${OUTFILEBASENAME}_${idx}.jld" $n $last
    let idx+=1
done
