"""
# Usage

    solve(FDDE::FDDEProblem, step_size, DelayPECE())

Using the delayed predictor-corrector method to solve the delayed fractional differential equation problem in the Caputo sense.

Capable of solving both single term FDDE and multiple FDDE, support time varying lags of course😋.

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
struct DelayPECE <: FDDEAlgorithm end
#FIXME: What if we have an FDDE with both variable order and time varying lag??
function solve(FDDE::FDDEProblem, step_size, ::DelayPECE)
    # If the delays are time varying, we need to specify single delay and multiple delay
    if  FDDE.constant_lags[1] isa Function
        # Here is the PECE solver for single time varying lag
        solve_fdde_with_single_lag(FDDE, step_size)
    elseif FDDE.constant_lags[1] isa AbstractArray{Function}
        # Here is the PECE solver for multiple time varying lags
        solve_fdde_with_multiple_lags(FDDE, step_size) #TODO: implement this
    # Varying order fractional delay differential equations
    elseif FDDE.order[1] isa Function
        if length(FDDE.constant_lags[1]) == 1
            # Here is the PECE solver for single lag with variable order
            solve_fdde_with_single_lag_and_variable_order(FDDE, step_size)
        else
            # Here is the PECE solver for multiple lags with variable order
            solve_fdde_with_multiple_lags_and_variable_order(FDDE, step_size)
        end
    else
        # If the delays are constant
        if length(FDDE.constant_lags[1]) == 1
            # Call the DelayPECE solver for single lag FDDE
            solve_fdde_with_single_lag(FDDE, step_size)
        else
            # Call the DelayPECE solver for multiple lags FDDE
            solve_fdde_with_multiple_lags(FDDE, step_size)
        end
    end
end

function solve_fdde_with_single_lag(FDDE::FDDEProblem, step_size)
    @unpack f, h, order, constant_lags, p, tspan = FDDE
    τ = constant_lags[1]
    iip = SciMLBase.isinplace(FDDE)
    T = tspan[2]
    t = collect(0:step_size:T)
    maxn = size(t, 1)
    yp = collect(Float64, 0:step_size:T+step_size)
    y = copy(t)
    y[1] = h(p, 0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            if iip
                tmp = zeros(length(yp[1]))
                f(tmp, y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
                yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, order, step_size)*tmp
            else
                yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, order, step_size)*f(y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
            end
        end
        yp[n+1] = yp[n+1]/gamma(order)+h(p, 0)

        y[n+1] = 0

        @fastmath @inbounds @simd for j=1:n
            if iip
                tmp = zeros(length(y[1]))
                f(tmp, y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
                y[n+1] = y[n+1]+a(j-1, n-1, order, step_size)*tmp
            else
                y[n+1] = y[n+1]+a(j-1, n-1, order, step_size)*f(y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
            end
        end

        if iip
            tmp = zeros(y[1])
            f(tmp, yp[n+1], v(h, n+1, τ, step_size, y, yp, t, p), p, t[n+1])
            y[n+1] = y[n+1]/gamma(order)+step_size^order*tmp/gamma(order+2)+h(p, 0)
        else
            y[n+1] = y[n+1]/gamma(order)+step_size^order*f(yp[n+1], v(h, n+1, τ, step_size, y, yp, t, p), p, t[n+1])/gamma(order+2)+h(p, 0)
        end
    end

    V = copy(t)
    @fastmath @inbounds @simd for n = 1:maxn-1
        # The delay term
        V[n] = v(h, n, τ, step_size, y, yp, t, p)
    end
    return V, y
end

function a(j::Int, n::Int, order, step_size)
    if j == n+1
        result = 1
    elseif j == 0
        result = n^(order+1)-(n-order)*(n+1)^order
    else
        result = (n-j+2)^(order+1) + (n-j)^(order+1) - 2*(n-j+1)^(order+1)
    end
    return result*step_size^order / (order*(order + 1))
end

#generalized_binomials(j, n, order::Number, step_size) = @. step_size^order/order*((n-j+1)^order - (n-j)^order)
generalized_binomials(j, n, order, step_size) = @. step_size^order/order*((n-j+1)^order - (n-j)^order)

function v(h, n, τ, step_size, y, yp, t, p)
    if typeof(τ) <: Function
        m = floor.(Int, τ.(t)/step_size)
        δ = m.-τ.(t)./step_size
        if τ(t[n]) >= (n-1)*step_size
            return h(p, (n-1)*step_size-τ(t[n]))
        else
            if m[n] > 1
                return δ[n]*y[n-m[n]+2] + (1-δ[n])*y[n-m[n]+1]
            elseif m[n] == 1
                return δ[n]*yp[n+1] + (1-δ[n])*y[n]
            end
        end
    else
        if τ >= (n-1)*step_size
            return h(p, (n-1)*step_size-τ)
        else
            m = floor(Int, τ/step_size)
            δ = m-τ/step_size

            if m>1
                return δ*y[n-m+2] + (1-δ)*y[n-m+1]
            elseif m == 1
                return δ*yp[n+1] + (1-δ)*y[n]
            end
        end
    end
end


function solve_fdde_with_multiple_lags(FDDE::FDDEProblem, step_size)
    @unpack f, h, order, constant_lags, p, tspan = FDDE
    τ = constant_lags[1]
    t = collect(0:step_size:tspan[2])
    maxn = length(t)
    yp = zeros(Float64, maxn)
    y = copy(t)
    y[1] = h(p, 0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, order, step_size)*f(t[j], y[j], multiv(h, j, τ, step_size, y, yp, p)...)
        end
        yp[n+1] = yp[n+1]/gamma(order)+h(p, 0)

        y[n+1] = 0

        for j=1:n
            y[n+1] = y[n+1]+multia(j-1, n-1, order, step_size)*f(t[j], y[j], multiv(h, j, τ, step_size, y, yp, p)...)
        end

        y[n+1] = y[n+1]/gamma(order)+step_size^order*f(t[n+1], yp[n+1], multiv(h, n+1, τ, step_size, y, yp, p)...)/gamma(order+2) + h(p, 0)
    end

    V = []
    for n = 1:maxn
        push!(V, multiv(h, n, τ, step_size, y, yp, p))
    end

    delayed = zeros(length(τ), length(V))
    for i=1:length(V)
        delayed[:, i] = V[i]
    end    
    return delayed, y
end

function multia(j::Int, n::Int, order, step_size)
    if j == n+1
        result = 1
    elseif j == 0
        result = n^(order+1)-(n-order)*(n+1)^order
    elseif j == n
        result = 2*(2^(order+1)-1)
    else
        result = (n-j+2)^(order+1) + (n-j)^(order+1) - 2*(n-j+1)^(order+1)
    end
    return result*step_size^order/(order*(order + 1))
end

function multiv(h, n, τ, step_size, y, yp, p)
    if maximum(τ) > n*step_size
        return h.(p, (n-1)*step_size.-τ)
    else
        m = floor.(Int, τ./step_size)
        δ = m.-τ./step_size

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

function solve_fdde_with_single_lag_and_variable_order(FDDE::FDDEProblem, step_size)
    @unpack f, order, h, constant_lags, p, tspan = FDDE
    iip = SciMLBase.isinplace(FDDE)
    order = order[1]
    τ = constant_lags[1]
    T = tspan[2]
    t = collect(0:step_size:T)
    maxn = size(t, 1)
    yp = collect(0:step_size:T+step_size)
    y = copy(t)
    y[1] = h(p, 0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            if iip
                tmp = zeros(length(yp[1]))
                f(tmp, y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
                yp[n+1] = yp[n+1] + generalized_binomials(j-1, n-1, order(t[n+1]), step_size)*tmp
            else
                yp[n+1] = yp[n+1] + generalized_binomials(j-1, n-1, order(t[n+1]), step_size)*f(y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
            end
        end
        yp[n+1] = yp[n+1]/gamma(order(t[n+1]))+h(p, 0)

        y[n+1] = 0

        @fastmath @inbounds @simd for j=1:n
            if iip
                tmp = zeros(length(y[1]))
                f(tmp, y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
                y[n+1] = y[n+1]+a(j-1, n-1, order(t[n+1]), step_size)*tmp
            else
                y[n+1] = y[n+1]+a(j-1, n-1, order(t[n+1]), step_size)*f(y[j], v(h, j, τ, step_size, y, yp, t, p), p, t[j])
            end
        end

        if iip
            tmp = zeros(length(y[1]))
            f(tmp, yp[n+1], v(h, n+1, τ, step_size, y, yp, t, p), p, t[n+1])
            y[n+1] = y[n+1]/gamma(order(t[n+1]))+step_size^order(t[n+1])*tmp/gamma(order(t[n+1])+2)+h(p, 0)
        else
            y[n+1] = y[n+1]/gamma(order(t[n+1]))+step_size^order(t[n+1])*f(yp[n+1], v(h, n+1, τ, step_size, y, yp, t, p), p, t[n+1])/gamma(order(t[n+1])+2)+h(p, 0)
        end
    end

    V = copy(t)
    @fastmath @inbounds @simd for n = 1:maxn-1
        V[n] = v(h, n, τ, step_size, y, yp, t, p)
    end
    return V, y
end


function solve_fdde_with_multiple_lags_and_variable_order(FDDE::FDDEProblem, step_size)
    @unpack f, h, order, constant_lags, p, tspan = FDDE
    τ = constant_lags[1]
    t = collect(0:step_size:tspan[2])
    maxn = length(t)
    yp = zeros(maxn)
    y = copy(t)
    y[1] = h(p, 0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+generalized_binomials(j-1, n-1, order(t[n+1]), step_size)*f(t[j], y[j], multiv(h, j, τ, step_size, y, yp, p)...)
        end
        yp[n+1] = yp[n+1]/gamma(order(t[n+1]))+h(p, 0)

        y[n+1] = 0

        for j=1:n
            y[n+1] = y[n+1]+multia(j-1, n-1, order(t[n+1]), step_size)*f(t[j], y[j], multiv(h, j, τ, step_size, y, yp, p)...)
        end

        y[n+1] = y[n+1]/gamma(order(t[n+1]))+step_size^order(t[n+1])*f(t[n+1], yp[n+1], multiv(h, n+1, τ, step_size, y, yp, p)...)/gamma(order(t[n+1])+2) + h(p, 0)
    end

    V = []
    for n = 1:maxn
        push!(V, multiv(h, n, τ, step_size, y, yp, p))
    end

    delayed = zeros(Float, length(τ), length(V))
    for i in eachindex(V)
        delayed[:, i] = V[i]
    end    
    return delayed, y
end