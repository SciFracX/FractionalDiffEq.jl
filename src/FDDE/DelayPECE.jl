"""
# Usage

    solve(FDDE::FDDEProblem, h, DelayPECE())

Using the delayed predictor-corrector method to solve the delayed fractional differential equation problem.

Capable of solving both single term FDDE and multiple FDDE.

### References

```tex
@article{Wang2013ANM,
  title={A Numerical Method for Delayed Fractional-Order Differential Equations},
  author={Zhen Wang},
  journal={J. Appl. Math.},
  year={2013},
  volume={2013},
  pages={256071:1-256071:7}
}
```
"""
struct DelayPECE <: FractionalDiffEqAlgorithm end

function solve(FDDE::FDDEProblem, h, ::DelayPECE)
    if length(FDDE.τ) > 1
        # Call the DelayPECE solver for multiple lags FDDE
        solve_fdde_with_multiple_lags(FDDE, h)
    else
        solve_fdde_with_single_lag(FDDE, h)
    end
end

function solve_fdde_with_single_lag(FDDE::FDDEProblem, h)
    @unpack f, ϕ, α, τ, tspan = FDDE
    T = tspan
    t = collect(0:h:T)
    maxn = size(t, 1)
    yp = collect(0:h:T+h)
    y = copy(t)
    y[1] = ϕ(0)

    @fastmath @inbounds @simd for n in 1:maxn-1
        yp[n+1] = 0
        @fastmath @inbounds @simd for j = 1:n
            yp[n+1] = yp[n+1]+b(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp))
        end
        yp[n+1] = yp[n+1]/gamma(α)+ϕ(0)

        y[n+1] = 0

        @fastmath @inbounds @simd for j=1:n
            y[n+1] = y[n+1]+a(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp))
        end

        y[n+1] = y[n+1]/gamma(α)+h^α*f(t[n+1], yp[n+1], v(ϕ, n+1, τ, h, y, yp))/gamma(α+2)+ϕ(0)
    end

    V = copy(t)
    @fastmath @inbounds @simd for n = 1:maxn-1
        V[n] = v(ϕ, n, τ, h, y, yp)
    end
    return V, y
end

function a(j, n, α, h)
    if j == n+1
        result = 1
    elseif j == 0
        result = n^(α+1)-(n-α)*(n+1)^α
    else
        result = (n-j+2)^(α+1) + (n-j)^(α+1) - 2*(n-j+1)^(α+1)
    end
    return result*h^α / (α*(α + 1))
end

function b(j, n, α, h)
    return h^α/α*((n-j+1)^α - (n-j)^α)
end

function v(ϕ, n, τ, h, y, yp)
    if typeof(τ) <: Function
        m = floor.(Int, τ/h)
        δ = @. m-τ/h
        if m[n] > 1
            return δ[n]*y[n-m[n]+2] + (1-δ[n])*y[n-m[n]+1]
        elseif m[n] == 1
            return δ[n]*yp[n+1] + (1-δ[n])*y[n]
        end
    else
        if τ >= (n-1)*h
            return ϕ((n-1)*h-τ)
        else
            m = floor(Int, τ/h)
            δ = m-τ/h

            if m>1
                return δ*y[n-m+2] + (1-δ)*y[n-m+1]
            elseif m == 1
                return δ*yp[n+1] + (1-δ)*y[n]
            end
        end
    end
end


function solve_fdde_with_multiple_lags(FDDE::FDDEProblem, h)
    @unpack f, ϕ, α, τ, tspan = FDDE
    t = collect(0:h:tspan)
    maxn = length(t)
    yp = zeros(maxn)
    y = copy(t)
    y[1] = ϕ(0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+b(j-1, n-1, α, h)*f(t[j], y[j], multiv(ϕ, j, τ, h, y, yp)...)
        end
        yp[n+1] = yp[n+1]/gamma(α)+ϕ(0)

        y[n+1] = 0

        for j=1:n
            y[n+1] = y[n+1]+multia(j-1, n-1, α, h)*f(t[j], y[j], multiv(ϕ, j, τ, h, y, yp)...)
        end

        y[n+1] = y[n+1]/gamma(α)+h^α*f(t[n+1], yp[n+1], multiv(ϕ, n+1, τ, h, y, yp)...)/gamma(α+2) + ϕ(0)
    end

    V = []
    for n = 1:maxn
        push!(V, multiv(ϕ, n, τ, h, y, yp))
    end

    delayed = zeros(length(τ), length(V))
    for i=1:length(V)
        delayed[:, i] = V[i]
    end    
    return delayed, y
end

function multia(j, n, α, h)
    if j == n+1
        result = 1
    elseif j == 0
        result = n^(α+1)-(n-α)*(n+1)^α
    elseif j == n
        result = 2*(2^(α+1)-1)
    else
        result = (n-j+2)^(α+1) + (n-j)^(α+1) - 2*(n-j+1)^(α+1)
    end
    return result*h^α/(α*(α + 1))
end

function multiv(ϕ, n, τ, h, y, yp)
    if maximum(τ) > n*h
        return ϕ.((n-1)*h.-τ)
    else
        m = floor.(Int, τ./h)
        δ = m.-τ./h

        function judge(m)
            temp1 = findall(x->x>1, m)
            temp2 = findall(x->x==1, m)#FIXME: Another case for x == 1

            result = zeros(length(m))
            if length(temp1) == length(m)
                for i=1:length(m)
                    result[i] = δ[i]*y[n-m[i]+2]+(1-δ[i])*y[n-m[i]+1]
                end
                return result
            end
        end
        return judge(m)
    end
end