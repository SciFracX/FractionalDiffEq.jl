"""
# Usage

    solve(FDDE::FDDEProblem, h, DelayPECE())

Using the delayed predictor-corrector method to solve the delayed fractional differential equation problem in the Caputo sense.

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

@inproceedings{Nagy2014NUMERICALSF,
  title={NUMERICAL SIMULATIONS FOR VARIABLE-ORDER FRACTIONAL NONLINEAR DELAY DIFFERENTIAL EQUATIONS},
  author={Abdelhameed M. Nagy and Taghreed Abdul Rahman Assiri},
  year={2014}
}

@inproceedings{Abdelmalek2019APM,
  title={A Predictor-Corrector Method for Fractional Delay-Differential System with Multiple Lags},
  author={Salem Abdelmalek and Redouane Douaifia},
  year={2019}
}
```
"""
struct DelayPECE <: AbstractFDEAlgorithm end
#FIXME: If we have a FDDE with both variable order and time varying lag??😂
function solve(FDDE::FDDEProblem, h, ::DelayPECE)
    # If the delays are time varying, we need to specify single delay and multiple delay
    if  typeof(FDDE.τ) <: Function
        # Here is the PECE solver for single time varying lag
        solve_fdde_with_single_lag(FDDE, h)
    elseif typeof(FDDE.τ) <: AbstractArray{Function}
        # Here is the PECE solver for multiple time varying lags
        solve_fdde_with_multiple_lags(FDDE, h) #TODO: implement this
    # Varying order fractional delay differential equations
    elseif typeof(FDDE.α) <: Function
        if length(FDDE.τ) == 1
            # Here is the PECE solver for single lag with variable order
            solve_fdde_with_single_lag_and_variable_order(FDDE, h)
        else
            # Here is the PECE solver for multiple lags with variable order
            solve_fdde_with_multiple_lags_and_variable_order(FDDE, h)
        end
    else
        # If the delays are constant
        if length(FDDE.τ) == 1
            # Call the DelayPECE solver for single lag FDDE
            solve_fdde_with_single_lag(FDDE, h)
        else
            # Call the DelayPECE solver for multiple lags FDDE
            solve_fdde_with_multiple_lags(FDDE, h)
        end
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

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp, t))
        end
        yp[n+1] = yp[n+1]/gamma(α)+ϕ(0)

        y[n+1] = 0

        @fastmath @inbounds @simd for j=1:n
            y[n+1] = y[n+1]+a(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp, t))
        end

        y[n+1] = y[n+1]/gamma(α)+h^α*f(t[n+1], yp[n+1], v(ϕ, n+1, τ, h, y, yp, t))/gamma(α+2)+ϕ(0)
    end

    V = copy(t)
    @fastmath @inbounds @simd for n = 1:maxn-1
        # The delay term
        V[n] = v(ϕ, n, τ, h, y, yp, t)
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

generalized_binomials(j, n, α, h) = h^α/α*((n-j+1)^α - (n-j)^α)

function v(ϕ, n, τ, h, y, yp, t)
    if typeof(τ) <: Function
        m = floor.(Int, τ.(t)/h)
        δ = m.-τ.(t)./h
        if τ(t[n]) >= (n-1)*h
            return ϕ((n-1)*h-τ(t[n]))
        else
            if m[n] > 1
                return δ[n]*y[n-m[n]+2] + (1-δ[n])*y[n-m[n]+1]
            elseif m[n] == 1
                return δ[n]*yp[n+1] + (1-δ[n])*y[n]
            end
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
            yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, α, h)*f(t[j], y[j], multiv(ϕ, j, τ, h, y, yp)...)
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

#########################For variable order FDDE###########################

function solve_fdde_with_single_lag_and_variable_order(FDDE::FDDEProblem, h)
    @unpack f, ϕ, α, τ, tspan = FDDE
    T = tspan
    t = collect(0:h:T)
    maxn = size(t, 1)
    yp = collect(0:h:T+h)
    y = copy(t)
    y[1] = ϕ(0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1] + generalized_binomials(j-1, n-1, α(t[n+1]), h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp, t))
        end
        yp[n+1] = yp[n+1]/gamma(α(t[n+1]))+ϕ(0)

        y[n+1] = 0

        @fastmath @inbounds @simd for j=1:n
            y[n+1] = y[n+1]+a(j-1, n-1, α(t[n+1]), h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp, t))
        end

        y[n+1] = y[n+1]/gamma(α(t[n+1]))+h^α(t[n+1])*f(t[n+1], yp[n+1], v(ϕ, n+1, τ, h, y, yp, t))/gamma(α(t[n+1])+2)+ϕ(0)
    end

    V = copy(t)
    @fastmath @inbounds @simd for n = 1:maxn-1
        V[n] = v(ϕ, n, τ, h, y, yp, t)
    end
    return V, y
end


function solve_fdde_with_multiple_lags_and_variable_order(FDDE::FDDEProblem, h)
    @unpack f, ϕ, α, τ, tspan = FDDE
    t = collect(0:h:tspan)
    maxn = length(t)
    yp = zeros(maxn)
    y = copy(t)
    y[1] = ϕ(0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, α(t[n+1]), h)*f(t[j], y[j], multiv(ϕ, j, τ, h, y, yp)...)
        end
        yp[n+1] = yp[n+1]/gamma(α(t[n+1]))+ϕ(0)

        y[n+1] = 0

        for j=1:n
            y[n+1] = y[n+1]+multia(j-1, n-1, α(t[n+1]), h)*f(t[j], y[j], multiv(ϕ, j, τ, h, y, yp)...)
        end

        y[n+1] = y[n+1]/gamma(α(t[n+1]))+h^α(t[n+1])*f(t[n+1], yp[n+1], multiv(ϕ, n+1, τ, h, y, yp)...)/gamma(α(t[n+1])+2) + ϕ(0)
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