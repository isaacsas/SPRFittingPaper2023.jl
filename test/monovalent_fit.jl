using SPRFitting, Test, LinearAlgebra
using Optimization
using OptimizationNLopt

# manufacture data
# parameters
kon = 1e-2
koff = 1e-1
CP = 8.5
toff = 150.0
tstop = 600.0
lb = [-8.0, -8.0, -8.0]
ub = [8.0, 8.0, 8.0]

method = NLopt.LD_LBFGS()
ad = Optimization.AutoForwardDiff()

tv = collect(range(0.0, tstop, length=601))
times = [tv, tv[1:2:end], tv[2:2:end]]
antibodyconcens = [25.0, 100.0, 150.0]
aligneddat = SPRFitting.AlignedData(times, totbndv, antibodyconcens, CP)

# forward solutions as data
u = (kon, koff, CP)
totbndv = similar(times)
for (i,abc) in enumerate(antibodyconcens)
    ts = times[i]
    totbnd = zeros(length(ts))
    p = (abc, toff)
    for (n,t) in enumerate(ts)
        totbnd[n] = SPRFitting.monovalent_total_bound(u, p, t)
    end

    totbndv[i] = totbnd
end

# using Plots
# p = plot()
# for (i,totbnd) in enumerate(totbndv)
#     plot!(p, times[i], totbndv[i], label=antibodyconcens[i])
# end
# plot(p)

# fit to forward solutions
u₀ = log10.(rand(3))
sol = monovalent_fit_spr_data(method, aligneddat, toff, u₀; ad, lb, ub, abstol=1e-12, reltol=1e-10)
ps = 10.0 .^ sol.u
@test norm(ps - [kon,koff,CP], Inf) < 1e-10