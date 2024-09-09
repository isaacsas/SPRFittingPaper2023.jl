using SPRFittingPaper2023, JLD

METADATAFILE = "/Users/isaacsas/data/2022-06-07 - FD11A_Data/surrogates/test_metadata.jld"
OUTFILEBASENAME = "/Users/isaacsas/data/2022-06-07 - FD11A_Data/surrogates/test_slice"
nfiles = 4   # number of files with the given basename
overwrite = true  # whether to overwrite an existing output file

# note this saves the data to a file named $OUTFILEBASENAME_merged.jld
merge_surrogate_slices(METADATAFILE, OUTFILEBASENAME, nfiles; force=overwrite)
