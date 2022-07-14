using Pkg

(length(ARGS) < 5) && error("Usage: julia pkg_path make_surrogate_parallel.jl metadatafile.jld outfile.jld startidx endidx")
pkgpath = ARGS[1]
metadatafile = ARGS[2]
isfile(metadatafile) || error("Error, metadata file named: $metadatafile note a valid file.")
outfile = ARGS[3]

# the range of parameter values to simulate
idxstart = parse(Int, ARGS[4])
idxend   = parse(Int, ARGS[5])

# BASEDIR = "/Users/isaacsas/data/2022-06-07 - FD11A_Data/surrogates"
# metadatafile = joinpath(BASEDIR, "test_metadata.jld")
# outfile = joinpath(BASEDIR, "test_slice_1.jld")
# idxstart = 1
# idxend = 3

########################## END INPUT #############################

Pkg.activate(pkgpath)

using SPRFitting, JLD

# save the surrogate slice
save_surrogate_slice(outfile, load(metadatafile), idxstart, idxend)