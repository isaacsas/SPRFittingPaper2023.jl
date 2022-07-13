using SPRFitting, JLD


(length(ARGS) < 4) && error("Usage: julia make_surrogate_parallel.jl metadatafile.jld outfile.jld startidx endidx")
metadatafile = ARGS[1]
isfile(metadatafile) || error("Error, metadata file named: $metadatafile note a valid file.")
outfile = ARGS[2]

# the range of parameter values to simulate
idxstart = parse(Int, ARGS[3])
idxend   = parse(Int, ARGS[4])

# BASEDIR = "/Users/isaacsas/data/2022-06-07 - FD11A_Data/surrogates"
# metadatafile = joinpath(BASEDIR, "test_metadata.jld")
# outfile = joinpath(BASEDIR, "test_slice_1.jld")
# idxstart = 1
# idxend = 3

########################## END INPUT #############################

surmetadata = load(metadatafile)

# build the slice of the surrogate
surrogate_data = build_surrogate_slice(surmetadata, idxstart, idxend)

# save the surrogate slice
save_surrogate_slice(outfile, surrogate_data, idxstart, idxend)