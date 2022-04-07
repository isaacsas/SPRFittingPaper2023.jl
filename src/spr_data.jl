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

# unit conversions 

# transforms RU to μM (is this really always right?)
RU_to_muM(x) = x/149/26795*1000000 

# transforms μM to (nm)⁻³ and vice versa
muM_to_inv_cubic_nm(x) = 6.023e-7 * x
inv_cubic_nm_to_muM(x) = x / 6.023e-7

# transforms RU to (nm)⁻³
RU_to_inv_cubic_nm(x) = x * (.6023/149/26795)