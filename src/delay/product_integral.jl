#=
# Usage

    solve(FDDE::FDDEProblem, dt, DelayPI())

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
=#
@concrete mutable struct DelayPIEXCache{iip, T}
    prob
    alg
    mesh
    u0
    order
    constant_lags
    p

    N
    y0
    y
    g
    b

    dt
    kwargs
end

Base.eltype(::DelayPIEXCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FDDEProblem, alg::DelayPIEX; dt=0.0, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    @unpack f, order, u0, h, tspan, p, constant_lags = prob
    τ = constant_lags[1]
    iip = SciMLBase.isinplace(prob)
    T = eltype(u0)
    t0 = tspan[1]; tfinal = tspan[2]
    N = ceil(Int, (tfinal-t0)/dt)
    mesh = t0 .+ dt*collect(0:N)

    order = order[1] # Only for commensurate order FDDE
    nn_al = collect(T, 0:N).^order[1]
    b = [0; nn_al[2:end].-nn_al[1:end-1]]/gamma(order+1)

    y0 = h(p, t0)
    length(y0) == 1 ? (y=[similar([y0]) for i=1:N+1]) : (y=[similar(y0) for i=1:N+1])

    length(y0) == 1 ? (g=[similar([y0]) for i=1:N+1]) : (g=[similar(y0) for i=1:N+1])
    length(y0) == 1 ? (y[1] = [y0]) : (y[1] = y0)

    return DelayPIEXCache{iip, T}(prob, alg, mesh, u0, order, τ, p, N, y0, y, g, b, dt, kwargs)
end

function SciMLBase.solve!(cache::DelayPIEXCache{iip, T}) where {iip, T}
    @unpack prob, alg, mesh, u0, order, constant_lags, p, N, y0, y, g, b, dt, kwargs = cache
    h_al = dt^order[1]
    τ = constant_lags
    l = length(u0)
    for n = 1:N
        tnm1 = mesh[n]
        if tnm1 <= τ
            y_nm1_tau = prob.h(p, tnm1-τ)
        else
            nm1_tau1 = floor(Int, (tnm1-τ)/dt)
            nm1_tau2 = ceil(Int, (tnm1-τ)/dt)
            if nm1_tau1 == nm1_tau2
                y_nm1_tau = y[nm1_tau1+1]
            else
                tt0 = mesh[nm1_tau1+1]
                tt1 = mesh[nm1_tau1+2]
                yy0 = y[nm1_tau1+1]
                yy1 = y[nm1_tau1+2]
                y_nm1_tau = ((tnm1-τ)-tt0)/(tt1-tt0)*yy1 + ((tnm1-τ)-tt1)/(tt0-tt1)*yy0
            end
        end

        if iip
            tmp = zeros(length(g[1]))
            prob.f(tmp, y[n], y_nm1_tau, p, tnm1)
            g[n] = tmp
        else
            length(u0) == 1 ? (tmp = prob.f(y[n], y_nm1_tau[1], p, tnm1)) : (tmp = prob.f(y[n], y_nm1_tau, p, tnm1))
            g[n] = tmp
        end
        f_mem = zeros(T, l)
        for j = 0:n-1
            @. f_mem = f_mem + g[j+1]*b[n-j+1]
        end
        @. y[n+1] = y0 + h_al*f_mem
    end

    return DiffEqBase.build_solution(prob, alg, mesh, y)
end