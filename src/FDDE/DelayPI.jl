"""
# Usage

    solve(FDDE::FDDEProblem, h, DelayPI())

Use explicit rectangular product integration algorithm to solve an FDDE problem.

### References

```tex
@article{2020,
   title={On initial conditions for fractional delay differential equations},
   ISSN={1007-5704},
   url={http://dx.doi.org/10.1016/j.cnsns.2020.105359},
   DOI={10.1016/j.cnsns.2020.105359},
   journal={Communications in Nonlinear Science and Numerical Simulation},
   author={Garrappa, Roberto and Kaslik, Eva},
   year={2020},
}
```
"""
struct DelayPI <: AbstractFDEAlgorithm end

function solve(FDDE::FDDEProblem, h, ::DelayPI)
    @unpack f, ϕ, α, τ, tspan = FDDE
    t0 = tspan[1]; T = tspan[2]
    N::Int = ceil(Int, (T-t0)/h)
    t = t0 .+ h*collect(0:N)

    nn_al = collect(Float64, 0:N).^α
    b = [0; nn_al[2:end].-nn_al[1:end-1]]/gamma(α+1)
    h_al = h^α

    y0 = ϕ(t0)
    y = zeros(Float64, N+1)

    g = zeros(Float64, N+1)
    y[1] = y0

    for n = 1:N
        tnm1 = t[n]
        if tnm1 <= τ
            y_nm1_tau = ϕ(tnm1-τ)
        else
            nm1_tau1 = floor(Int, (tnm1-τ)/h)
            nm1_tau2 = ceil(Int, (tnm1-τ)/h)
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

        g[n] = f(tnm1, y[n], y_nm1_tau)
        f_mem = 0
        for j = 0:n-1
            f_mem = f_mem + g[j+1]*b[n-j+1]
        end
        y[n+1] = y0 + h_al*f_mem
    end
    return y
end