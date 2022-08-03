"""
$(TYPEDEF)

The aligned data from a set of SPR experiments varying antibody concentrations.

# Fields
$(FIELDS)
"""
Base.@kwdef struct AlignedData
    """times[i] gives the SPR data timepoints for antibodyconcen[i]"""
    times::Vector{Vector{Float64}}
    """refdata[i] gives the SPR simulation data for antibodyconcens[i] """
    refdata::Vector{Vector{Float64}}
    """Antibody concentrations for the SPR data set in μM"""
    antibodyconcens::Vector{Float64}
    """Antigen concentration for the SPR data set in μM"""
    antigenconcen::Float64
end

"""
    get_aligned_data(fname, antigenconcen)

Read the given CSV file representing aligned SPR data. Returns an `AlignedData` structure
storing the given SPR experiment data.

Notes:
- `antigenconcen = nothing` the concentration of antigen used in the experiments in units of
  μM. If not set assumes the CSV filename has the form ""text-AGCCONCEN_aligned.csv"", and
  parses the number antigen concentration from AGCCONCEN.
"""
function get_aligned_data(fname, antigenconcen=nothing)

    # assumes fname = "text-AGCCONCEN_aligned.csv"
    if antigenconcen === nothing
        agcstr = split(split(fname, '-')[end], '_')[end-1]
        antigenconcen = parse(Float64, agcstr)
    end

    cf = CSV.File(fname)
    headers = Tables.columnnames(cf)
    nabcs   = div(length(headers), 2)
    hstrs = map(first ∘ Base.Fix2(split, "_") ∘ String, headers)
    antibodyconcens = parse.(Float64, hstrs[2:2:end])

    cols = Tables.columns(cf)
    idx  = 1
    times = Vector{Vector{Float64}}(undef, nabcs)
    refdata = Vector{Vector{Float64}}(undef, nabcs)
    for (i,colname) in enumerate(headers)
        col = Tables.getcolumn(cols, colname)
        if isodd(i)
            times[idx] = collect(skipmissing(col))
        else
            refdata[idx] = collect(skipmissing(col))
            (length(refdata[idx]) == length(times[idx])) || error("Found a column of SPR data times with a different length than the assoicated column of SPR curve data.")
            idx += 1
        end

    end

    return AlignedData(times, refdata, antibodyconcens, antigenconcen)
end
