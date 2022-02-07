"""
```tex
@article{2020,
   title={On initial conditions for fractional delay differential equations},
   volume={90},
   ISSN={1007-5704},
   url={http://dx.doi.org/10.1016/j.cnsns.2020.105359},
   DOI={10.1016/j.cnsns.2020.105359},
   journal={Communications in Nonlinear Science and Numerical Simulation},
   publisher={Elsevier BV},
   author={Garrappa, Roberto and Kaslik, Eva},
   year={2020},
   month={Nov},
   pages={105359}
}
```
"""
struct DelayPI <: FractionalDiffEqAlgorithm end

function solve(FDDE::FDDEProblem, T, h, ::DelayPI)
    g, ϕ, α, τ, t0 = FDDE.f, FDDE.ϕ, FDDE.α, FDDE.τ, FDDE.t0
    N = Int64(ceil((T-t0)/h))
    t = t0 .+ h*collect(0:1:N)

    nn_al = collect(0:1:N).^α
    b = [0; nn_al[2:end].-nn_al[1:end-1]]/gamma(α+1)
    h_al = h^α

    y0 = ϕ(t0)
    y = zeros(N+1)

    f = zeros(N+1)
    y[1] = y0

    for n = 1:N
        tnm1 = t[n]
        if tnm1 <= τ
            y_nm1_tau = ϕ(tnm1-τ)
        else
            nm1_tau1 = Int64(floor((tnm1-τ)/h))
            nm1_tau2 = Int64(ceil((tnm1-τ)/h))
            if nm1_tau1 == nm1_tau2
                y_nm1_tau = y[nm1_tau1+1]
            else
                tt0 = t[nm1_tau1+1]
                tt1 = t[nm1_tau1+2]
                yy0 = y[nm1_tau1+1]
                yy1 = y[nm1_tau1+2]
                y_nm1_tau = ((tnm1-τ)-tt0)/(tt1-tt0)*yy1 + ((tnm1-τ)-tt1)/(tt0-tt1)*yy0
            end
        end

        f[n] = g(tnm1, y[n], y_nm1_tau)
        f_mem = 0
        for j = 0:n-1
            f_mem = f_mem + f[j+1]*b[n-j+1]
        end
        y[n+1] = y0 + h_al*f_mem
    end
    return y
end