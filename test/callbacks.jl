using SPRFitting, Statistics, StableRNGs, Test
rng = StableRNG(12345)

# check calculating mean and variance ok
nsims = 10
x   = rand(rng,3,nsims)
tbo = TotalBoundOutputter(1)
biopars = (CP = 1.0,)
numpars = (N=1, nsims=size(x,2))
for xcol in eachcol(x)
    tbo(0.0, xcol, biopars, numpars)    
end
tbo(biopars,numpars)
m   = tbo.bindcnt[1]
sse = tbo.bindcntsse[1]
xvec = x[2,:] + x[3,:]
@test mean(xvec) ≈ m
@test (std(xvec) / sqrt(nsims) / mean(xvec)) ≈ sse