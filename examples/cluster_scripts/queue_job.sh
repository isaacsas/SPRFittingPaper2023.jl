#! /bin/bash

PKGPATH="/project/fpkmc3d/.julia_sai/dev/SPRFitting"
METADATAFILE="/project/fpkmc3d/surrogates/tests/test_metadata.jld"
OUTFILEBASENAME="/project/fpkmc3d/surrogates/tests/test_slice"
JULIASCRIPT="${PKGPATH}/examples/make_surrogate_parallel.jl"
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