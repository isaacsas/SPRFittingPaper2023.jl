@doc raw"""
    monovalent_total_bound(u,p,t)

Evaluates the monovalent model solution at time `t`.

The amount of bound antigen at time `t` is
```math
\begin{aligned}
A(t) = \begin{cases}
\frac{k_{\text{on}} A_b}{k_{\text{on}} A_b + k_{\text{off}}} * \left(1 - e^{-(k_{\text{on}} A_b + k_{\text{off}}) t} \right), & t \leq t_{\text{off}} \\
A(t_{\text{off}}) e^{-k_{\text{off}} (t - t_{\text{off}})}, & t > t_{\text{off}}
\end{cases}
\end{aligned}
```
where ``A_b`` is the amount of antibody and ``t_{\text{off}}`` is the time antibody stops
being injected.

Notes:
- `u = [kon, koff, CP]`
- `p = [Antibody, toff]`, where `toff` is the time at which antibody is removed.
"""
function monovalent_total_bound(u,p,t)
    kon, koff, CP = u
    antibody, toff = p
    δ = (kon*antibody + koff)
    A = (t <= toff) ? (1 - exp(-δ*t)) : (1 - exp(-δ*toff)) * exp(-koff * (t-toff))
    return (CP * kon * antibody * A) / δ
end

"""
    monovalent_sprdata_error(optpars, aligned_data::AlignedData, toff)

Evaluates the ℓ₂ error between the SPR data and the monovalent model for the given
optimization parameters.

Notes:
- `optpars = [log₁₀(kon), log₁₀(koff), log₁₀(CP)]`
- `toff` should be the time that antibodies were removed from the system.
"""
function monovalent_sprdata_error(optpars, pars)
    @unpack times, refdata, antibodyconcens = pars.aligneddat
    kon = 10.0 ^ optpars[1]
    koff = 10.0 ^ optpars[2]
    CP = 10.0 ^ optpars[3]
    err = 0.0
    @inbounds for (j,abc) in enumerate(antibodyconcens)
        u = (kon, koff, CP)
        p = (abc, pars.toff)
        sprdata = refdata[j]
        @inbounds for (i,t) in enumerate(times[j])
            newerr = monovalent_total_bound(u, p, t) - sprdata[i]
            err += newerr * newerr
        end
    end

    # return the ℓ₂ error
    return err
end

function monovalent_fit_spr_data(optimiser::Optimiser, aligneddat::AlignedData, toff, u₀)
    @unpack method, ad, lb, ub, probkwargs, solverkwargs = optimiser

    optfun = if ad === nothing
        OptimizationFunction(monovalent_sprdata_error)
    else
        OptimizationFunction(monovalent_sprdata_error, ad)
    end

    pars = (; aligneddat, toff)
    prob = OptimizationProblem(optfun, u₀, pars; lb, ub, probkwargs...)
    solve(prob, method; solverkwargs...)
end