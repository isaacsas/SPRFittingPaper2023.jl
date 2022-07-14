#! /bin/bash

PKGPATH="/project/fpkmc3d/.julia_sai/dev/SPRFitting"
METADATAFILE="/Users/isaacsas/data/2022-06-07 - FD11A_Data/surrogates/test_metadata.jld"
OUTFILEBASENAME="/Users/isaacsas/data/2022-06-07 - FD11A_Data/surrogates/test_slice"
JULIASCRIPT="make_surrogate_parallel.jl"
NUMPARAMSETS=8
NCPUS=4

let IDXINCREMENT=$NUMPARAMSETS/$NCPUS
let dn=$IDXINCREMENT-1
let idx=1
for (( n=1;n<=$NUMPARAMSETS;n+=$IDXINCREMENT)); do
    let last=$n+$dn
    qsub run_sim.sh "$JULIASCRIPT" "$PKGPATH" "$METADATAFILE" "${OUTFILEBASENAME}_${idx}.jld" $n $last
    let idx+=1
done