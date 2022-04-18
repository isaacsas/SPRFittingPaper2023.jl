"""
$(TYPEDEF)

The aligned data from a set of SPR experiments varying antibody concentrations.

# Fields
$(FIELDS)
"""
Base.@kwdef struct AlignedData
    """Times SPR data is defined at"""
    times::Vector{Float64}
    """length(times) by length(antibodyconcens) matrix of SPR data"""
    refdata::Matrix{Float64}
    """refdata with each nan replaced by zero"""
    refdata_nonan::Matrix{Float64}
    """Antibody concentrations for the SPR data set in μM"""
    antibodyconcens::Vector{Float64}
    """Antigen concentration for the SPR data set in μM"""
    antigenconcen::Float64
end

"""
    get_aligned_data(fname, antigenconcen)

Read the given CSV file representing aligned SPR data. Returns an `AlignedData`
structure storing the given SPR experiment data.

Notes:
- `antigenconcen = ` the concentration of antigen used in the experiments in
  units of μM.
"""
function get_aligned_data(fname, antigenconcen)

    data, header = readdlm(fname,',', header=true)
    
    antibodyconcens = parse.(Float64, header[2:end])
    times   = data[:,1]    
    refdata = data[:,2:end]

    # replace nan's by zeros
    refdata_nonan = map(x -> isnan(x) ? 0.0 : x, refdata)
    
    return AlignedData(times, refdata, refdata_nonan, antibodyconcens, antigenconcen)
end
