using SPRFitting, JLD

# the file in which you saved the metadata
METADATAFILE = "/project/fpkmc3d/surrogates/rebuild_lut_highandlow/surrogate_metadata.jld"

# the basename of the files with the output from each job
OUTFILEBASENAME = "/project/fpkmc3d/surrogates/rebuild_lut_highandlow/surrogate_slice"

nfiles = 900   # number of files with the given basename
overwrite = false  # whether to overwrite an existing output file

# note this saves the data to a file named $OUTFILEBASENAME_merged.jld
merge_surrogate_slices(METADATAFILE, OUTFILEBASENAME, nfiles; force=overwrite)
