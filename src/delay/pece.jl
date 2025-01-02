@concrete mutable struct DelayPECECache{iip, T}
    prob
    alg
    mesh
    u0
    order
    constant_lags
    p
    y
    yp
    L

    dt
    kwargs
end
#=
#FIXME: What if we have an FDDE with both variable order and time varying lag??
function solve(FDDE::FDDEProblem, dt, ::DelayPECE)
    # If the delays are time varying, we need to specify single delay and multiple delay
    if  FDDE.constant_lags[1] isa Function
        # Here is the PECE solver for single time varying lag
        solve_fdde_with_single_lag(FDDE, dt)
    elseif FDDE.constant_lags[1] isa AbstractArray{Function}
        # Here is the PECE solver for multiple time varying lags
        solve_fdde_with_multiple_lags(FDDE, dt) #TODO: implement this
    # Varying order fractional delay differential equations
    elseif FDDE.order[1] isa Function
        if length(FDDE.constant_lags[1]) == 1
            # Here is the PECE solver for single lag with variable order
            solve_fdde_with_single_lag_and_variable_order(FDDE, dt)
        else
            # Here is the PECE solver for multiple lags with variable order
            solve_fdde_with_multiple_lags_and_variable_order(FDDE, dt)
        end
    else
        # If the delays are constant
        if length(FDDE.constant_lags[1]) == 1
            # Call the DelayPECE solver for single lag FDDE
            solve_fdde_with_single_lag(FDDE, dt)
        else
            # Call the DelayPECE solver for multiple lags FDDE
            solve_fdde_with_multiple_lags(FDDE, dt)
        end
    end
end
=#

function SciMLBase.__init(prob::FDDEProblem, alg::DelayPECE; dt = 0.0, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    (; f, h, order, u0, constant_lags, p, tspan) = prob
    τ = constant_lags[1]
    iip = SciMLBase.isinplace(prob)
    t0 = tspan[1]
    tfinal = tspan[2]
    T = eltype(u0)
    N = Int(cld(tfinal - t0, dt))
    L = length(collect(-τ:dt:t0))
    mesh = collect(range(t0; stop = tfinal, length = N + 1))
    yp = _generate_similar_array(u0, N+1, u0)
    y = _generate_similar_array(u0, L+N+1, u0)
    y[L:-1:1] .= ifelse((length(u0) == 1), map(x->[h(p, x)], collect(-τ:dt:t0)), map(x->h(p, x), collect(-τ:dt:(t0))))
    y[L+1] = ifelse((length(u0) == 1), [u0], u0)#_generate_similar_array(u0, 1, h(p, 0))

    return DelayPECECache{iip, T}(prob, alg, mesh, u0, order, τ, p, y, yp, L, dt, kwargs)
end

@inline _generate_similar_array(u0, N, src) = ifelse(length(u0) == 1, [[zero(u0)] for _ in 1:N], [zero(src) for _ in 1:N])

@inline function OrderWrapper(order, t)
    if order isa Function
        return order(t)
    else
        return order
    end
end

function SciMLBase.solve!(cache::DelayPECECache{iip, T}) where {iip, T}
    if length(cache.constant_lags) == 1
        return solve_fdde_with_single_lag(cache)
    else
        return solve_fdde_with_multiple_lags(cache)
    end
end

function solve_fdde_with_single_lag(cache::DelayPECECache{iip, T}) where {iip, T}
    (; prob, alg, mesh, u0, order, L, p, dt) = cache
    N = length(mesh)
    l = length(u0)
    initial = ifelse(l == 1, [u0], u0)
    tmp = zeros(T, l)
    for i in 0:N-2
        order = OrderWrapper(order, mesh[i + 1])
        # compute the yp part
        for j in 0:i
            if iip
                prob.f(tmp, cache.y[L+j+1], v(cache, j), p, mesh[j + 1])
                cache.yp[i + 1] = cache.yp[i + 1] +
                                  generalized_binomials(j, i, order, dt) * tmp
            else
                length(u0) == 1 ?
                (tmp = prob.f(cache.y[L+j+1][1], v(cache, j)[1], p, mesh[j + 1])) :
                (tmp = prob.f(cache.y[L+j+1], v(cache, j), p, mesh[j+1]))
                @. cache.yp[i + 1] = cache.yp[i + 1] +
                                     generalized_binomials(j, i, order, dt) * tmp
            end
        end
        cache.yp[i + 1] = cache.yp[i + 1] / gamma(order) + initial

        # compute the a part
        for j in 0:i
            if iip
                prob.f(tmp, cache.y[L + j + 1], v(cache, j), p, mesh[j+1])
                cache.y[L + i + 2] += a(j, i, order, dt) * tmp
            else
                length(u0) == 1 ?
                (tmp = prob.f(cache.y[L + j + 1][1], v(cache, j)[1], p, mesh[j+1])) :
                (tmp = prob.f(cache.y[L + j + 1], v(cache, j), p, mesh[j]))
                @. cache.y[L + i + 2] += a(j, i, order, dt) * tmp
            end
        end

        if iip
            prob.f(tmp, cache.yp[i + 1], v(cache, i + 1), p, mesh[i + 2])
            cache.y[L + i + 2] = cache.y[L + i + 2] / gamma(order) +
                             dt^order * tmp / gamma(order + 2) + initial
        else
            length(u0) == 1 ?
            (tmp = prob.f(cache.yp[i + 1][1], v(cache, i + 1)[1], p, mesh[i + 2])) :
            (tmp = prob.f(cache.yp[i + 1], v(cache, i + 1), p, mesh[i + 2]))
            @. cache.y[L + i + 2] = cache.y[L + i + 2] / gamma(order) +
                                dt^order * tmp / gamma(order + 2) + initial
        end
    end

    u = cache.y[L+1:end]

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end

function a(j::I, n::I, order, dt) where {I <: Integer}
    if j == n
        result = 1
    elseif j == 0
        result = n^(order + 1) - (n - order) * (n + 1)^order
    else
        result = (n - j + 2)^(order + 1) + (n - j)^(order + 1) - 2 * (n - j + 1)^(order + 1)
    end
    return result * dt^order / (order * (order + 1))
end

function generalized_binomials(j, n, order, dt)
    @. dt^order / order * ((n - j + 1)^order - (n - j )^order)
end

function v(cache::DelayPECECache{iip, T}, n) where {iip, T}
    (; prob, mesh, dt, constant_lags, L, p) = cache
    τ = constant_lags
    if typeof(τ) <: Function
        m = floor.(Int, τ.(mesh) / dt)
        δ = m .- τ.(mesh) ./ dt
        if τ(mesh[n]) >= (n - 1) * dt
            return prob.h(p, (n - 1) * dt - τ(mesh[n]))
        else
            if m[n] > 1
                return δ[n] * cache.y[n - m[n] + 2] + (1 - δ[n]) * cache.y[n - m[n] + 1]
            elseif m[n] == 1
                return δ[n] * cache.yp[n + 1] + (1 - δ[n]) * cache.y[n]
            end
        end
    else
        if isinteger(τ/dt)#y-τ fall in grid point
            # case when δ == 0
            if τ/dt < n
                return cache.y[L + n - Int(τ/dt) + 1]
            else
                if length(prob.u0) == 1
                    return [prob.h(p, n * dt - τ)]
                else
                    return prob.h(p, n * dt - τ)
                end
            end
        else
            m = floor(Int, τ / dt)
            δ = m - τ / dt
            # interpolate by the two nearest points
            if m > 1
                return δ * cache.y[n - m + L + 2] + (1 - δ) * cache.y[n - m + L + 1]
            elseif m == 1 # h is larger than τ
                return δ * cache.yp[n + 1] + (1 - δ) * cache.y[L + n + 1]
            end
        end
    end
end

function solve_fdde_with_multiple_lags(cache::DelayPECECache)
    (; f, h, order, constant_lags, p, tspan) = cache
    τ = constant_lags[1]
    mesh = collect(0:dt:tspan[2])
    maxn = length(mesh)
    yp = zeros(Float64, maxn)
    y = copy(mesh)
    y[1] = h(p, 0)

    for n in 1:(maxn - 1)
        yp[n + 1] = 0
        for j in 1:n
            yp[n + 1] = yp[n + 1] +
                        generalized_binomials(j - 1, n - 1, order, dt) *
                        f(mesh[j], y[j], multiv(h, j, τ, dt, y, yp, p)...)
        end
        yp[n + 1] = yp[n + 1] / gamma(order) + h(p, 0)

        y[n + 1] = 0

        for j in 1:n
            y[n + 1] = y[n + 1] +
                       multia(j - 1, n - 1, order, dt) *
                       f(mesh[j], y[j], multiv(h, j, τ, dt, y, yp, p)...)
        end

        y[n + 1] = y[n + 1] / gamma(order) +
                   dt^order *
                   f(mesh[n + 1], yp[n + 1], multiv(h, n + 1, τ, dt, y, yp, p)...) /
                   gamma(order + 2) +
                   h(p, 0)
    end

    V = []
    for n in 1:maxn
        push!(V, multiv(h, n, τ, dt, y, yp, p))
    end

    delayed = zeros(length(τ), length(V))
    for i in 1:length(V)
        delayed[:, i] = V[i]
    end
    return delayed, y
end

function multia(j::Int, n::Int, order, dt)
    if j == n + 1
        result = 1
    elseif j == 0
        result = n^(order + 1) - (n - order) * (n + 1)^order
    elseif j == n
        result = 2 * (2^(order + 1) - 1)
    else
        result = (n - j + 2)^(order + 1) + (n - j)^(order + 1) - 2 * (n - j + 1)^(order + 1)
    end
    return result * dt^order / (order * (order + 1))
end

function multiv(h, n, τ, dt, y, yp, p)
    if maximum(τ) > n * dt
        return h.(p, (n - 1) * dt .- τ)
    else
        m = floor.(Int, τ ./ dt)
        δ = m .- τ ./ dt

        function judge(m)
            temp1 = findall(x -> x > 1, m)
            temp2 = findall(x -> x == 1, m)#FIXME: Another case for x == 1

            result = zeros(length(m))
            if length(temp1) == length(m)
                for i in 1:length(m)
                    result[i] = δ[i] * y[n - m[i] + 2] + (1 - δ[i]) * y[n - m[i] + 1]
                end
                return result
            end
        end
        return judge(m)
    end
end

#########################For variable order FDDE###########################

function solve_fdde_with_single_lag_and_variable_order(FDDE::FDDEProblem, dt)
    (; f, order, h, constant_lags, p, tspan) = FDDE
    iip = SciMLBase.isinplace(FDDE)
    order = order[1]
    τ = constant_lags[1]
    tfinal = tspan[2]
    mesh = collect(0:dt:tfinal)
    maxn = size(mesh, 1)
    yp = collect(0:dt:(tfinal + dt))
    y = copy(mesh)
    y[1] = h(p, 0)

    for n in 1:(maxn - 1)
        yp[n + 1] = 0
        for j in 1:n
            if iip
                tmp = zeros(length(yp[1]))
                f(tmp, y[j], v(h, j, τ, dt, y, yp, mesh, p), p, mesh[j])
                yp[n + 1] = yp[n + 1] +
                            generalized_binomials(j - 1, n - 1, order(mesh[n + 1]), dt) *
                            tmp
            else
                yp[n + 1] = yp[n + 1] +
                            generalized_binomials(j - 1, n - 1, order(mesh[n + 1]), dt) *
                            f(y[j], v(h, j, τ, dt, y, yp, mesh, p), p, mesh[j])
            end
        end
        yp[n + 1] = yp[n + 1] / gamma(order(mesh[n + 1])) + h(p, 0)

        y[n + 1] = 0

        @fastmath @inbounds @simd for j in 1:n
            if iip
                tmp = zeros(length(y[1]))
                f(tmp, y[j], v(h, j, τ, dt, y, yp, mesh, p), p, mesh[j])
                y[n + 1] = y[n + 1] + a(j - 1, n - 1, order(mesh[n + 1]), dt) * tmp
            else
                y[n + 1] = y[n + 1] +
                           a(j - 1, n - 1, order(mesh[n + 1]), dt) *
                           f(y[j], v(h, j, τ, dt, y, yp, mesh, p), p, mesh[j])
            end
        end

        if iip
            tmp = zeros(length(y[1]))
            f(tmp, yp[n + 1], v(h, n + 1, τ, dt, y, yp, mesh, p), p, mesh[n + 1])
            y[n + 1] = y[n + 1] / gamma(order(mesh[n + 1])) +
                       dt^order(mesh[n + 1]) * tmp / gamma(order(mesh[n + 1]) + 2) +
                       h(p, 0)
        else
            y[n + 1] = y[n + 1] / gamma(order(mesh[n + 1])) +
                       dt^order(mesh[n + 1]) *
                       f(yp[n + 1], v(h, n + 1, τ, dt, y, yp, mesh, p), p, mesh[n + 1]) /
                       gamma(order(mesh[n + 1]) + 2) +
                       h(p, 0)
        end
    end

    V = copy(mesh)
    @fastmath @inbounds @simd for n in 1:(maxn - 1)
        V[n] = v(h, n, τ, dt, y, yp, mesh, p)
    end
    return V, y
end

function solve_fdde_with_multiple_lags_and_variable_order(FDDE::FDDEProblem, dt)
    (; f, h, order, constant_lags, p, tspan) = FDDE
    τ = constant_lags[1]
    mesh = collect(0:dt:tspan[2])
    maxn = length(mesh)
    yp = zeros(maxn)
    y = copy(mesh)
    y[1] = h(p, 0)

    for n in 1:(maxn - 1)
        yp[n + 1] = 0
        for j in 1:n
            yp[n + 1] = yp[n + 1] +
                        generalized_binomials(j - 1, n - 1, order(mesh[n + 1]), dt) *
                        f(mesh[j], y[j], multiv(h, j, τ, dt, y, yp, p)...)
        end
        yp[n + 1] = yp[n + 1] / gamma(order(mesh[n + 1])) + h(p, 0)

        y[n + 1] = 0

        for j in 1:n
            y[n + 1] = y[n + 1] +
                       multia(j - 1, n - 1, order(mesh[n + 1]), dt) *
                       f(mesh[j], y[j], multiv(h, j, τ, dt, y, yp, p)...)
        end

        y[n + 1] = y[n + 1] / gamma(order(mesh[n + 1])) +
                   dt^order(mesh[n + 1]) *
                   f(mesh[n + 1], yp[n + 1], multiv(h, n + 1, τ, dt, y, yp, p)...) /
                   gamma(order(mesh[n + 1]) + 2) +
                   h(p, 0)
    end

    V = []
    for n in 1:maxn
        push!(V, multiv(h, n, τ, dt, y, yp, p))
    end

    delayed = zeros(Float, length(τ), length(V))
    for i in eachindex(V)
        delayed[:, i] = V[i]
    end
    return delayed, y
end
