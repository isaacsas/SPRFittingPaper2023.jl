using SPRFittingPaper2023, JLD

METADATAFILE = "/project/fpkmc3d/surrogates/highandlow_widerkoff/surrogate_metadata.jld"
OUTFILEBASENAME = "/project/fpkmc3d/surrogates/highandlow_widerkoff/surrogate_slice"
nfiles = 1200   # number of files with the given basename
overwrite = false  # whether to overwrite an existing output file

# note this saves the data to a file named $OUTFILEBASENAME_merged.jld
merge_surrogate_slices(METADATAFILE, OUTFILEBASENAME, nfiles; force=overwrite)
