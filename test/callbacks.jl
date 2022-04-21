using SPRFitting, Statistics, StableRNGs, OnlineStats, Test
using SPRFitting: means, sems
rng = StableRNG(12345)

# check calculating mean and variance ok
nsims = 10
x   = rand(rng,3,nsims)
tbo = TotalBoundOutputter(1)
biopars = (CP = 1.0,)
simpars = (N=1, nsims=size(x,2))
for xcol in eachcol(x)
    tbo(0.0, xcol, biopars, simpars)    
end
tbo(biopars, simpars)
m   = means(tbo)[1]
sse = sems(tbo)[1] / means(tbo)[1]
xvec = x[2,:] + x[3,:]
@test mean(xvec) ≈ m
@test (std(xvec) / sqrt(nsims) / mean(xvec)) ≈ sse

# test variance terminator
biopars = BioPhysParams(; kon = .001, koff=.01, konb=.01, reach=10.0, 
                          antigenconcen=SPRFitting.DEFAULT_SIM_ANTIGENCONCEN)
simpars = SimParams(; tstop_AtoB=150.0, dt=1.0)
tbo     = TotalBoundOutputter(length(simpars.tsave))
vt      = SPRFitting.VarianceTerminator(; ssetol=.025, maxsims=4000)
run_spr_sim!(tbo, biopars, simpars, vt)
function ssem(tbo)
    sqnsims = sqrt(nobs(tbo.bindcnt[1]))
    v = sems(tbo) ./ means(tbo)
    for i in eachindex(v)
        isnan(v[i]) && (v[i] = 0.0)
    end
    v
end
@show maximum(ssem(tbo)), vt.ssetol, nobs(tbo.bindcnt[1])
@test all(std(o)/sqrt(nobs(o)) <= vt.ssetol*mean(o) for o in tbo.bindcnt) || (nobs(tbo.bindcnt[1]) == vt.maxsims)
