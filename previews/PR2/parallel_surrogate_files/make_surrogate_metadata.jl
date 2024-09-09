using SPRFitting, JLD

# directory these scripts and generated surrogate data will be in 
BASEDIR = "/project/fpkmc3d/surrogates/rebuild_lut_highandlow"
mkpath(BASEDIR)

# within the directory the file to save the metadata within
metadatafile = joinpath(BASEDIR, "surrogate_metadata.jld")

# first we set biophysical parameters for the forward simulations
# that build the surrogate.
logkon_range  = (-5.0, 2.0)
logkoff_range = (-4.0, -1.0)
logkonb_range = (-3.0, 1.5)
reach_range   = (2.0, 35.0)

# now we set the numerical parameters, see SimParams for others
tstop      = 600.0               # simulation end time
tsavelen   = 601                 # number of times to save 
tstop_AtoB = 150.0               # time to turn off the surrogate

# size of the surrogate in each coordinate: nkon,nkoff,nkonb,nreach,tsavelen
surrogate_size = (42,30,30,30,tsavelen)

########################## END INPUT #############################

# make the output directory if not already present
DIR = dirname(metadatafile)
mkpath(DIR)

surpars = SurrogateParams(; logkon_range, logkoff_range, logkonb_range, reach_range)
tsave   = collect(range(0.0, tstop, tsavelen))
simpars = SimParams(; tstop, tstop_AtoB, tsave)

# save the surrogate
save_surrogate_metadata(metadatafile, surrogate_size, surpars, simpars)