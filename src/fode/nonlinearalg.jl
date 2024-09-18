@concrete mutable struct NonlinearAlgCache{iip, T}
    prob
    alg
    mesh
    u0
    order
    halpha
    p
    N
    problem_size

    gen
    W
    z
    y
    L0
    memory_effect
    kwargs
end

Base.eltype(::NonlinearAlgCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(
        prob::FODEProblem, alg::NonLinearAlg; dt = 0.0, L0 = 1e10, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    prob = _is_need_convert!(prob)
    (; f, order, u0, tspan, p) = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    T = eltype(u0)
    iip = isinplace(prob)
    mesh = collect(t0:dt:tfinal)
    problem_size = length(u0)
    N = round(Int, (tfinal - t0) / dt) + 1
    gen = genfun(1)
    gen = gen[:]
    u0 = u0[:]
    halpha = dt .^ order
    z = zeros(T, problem_size, N)

    # All of the min(N, L0+1) is to set the memory effect.
    memory_effect = Int(min(N, L0 + 1))
    W = zeros(T, problem_size, memory_effect) #Initializing W a problem_size*N matrix
    return NonlinearAlgCache{iip, T}(prob, alg, mesh, u0, order, halpha, p, N, problem_size,
        gen, W, z, similar(z), L0, memory_effect, kwargs)
end

function SciMLBase.solve!(cache::NonlinearAlgCache{iip, T}) where {iip, T}
    (; prob, alg, mesh, u0, order, halpha, p, N, problem_size, gen, W, z, y, L0, memory_effect, kwargs) = cache

    @fastmath @inbounds @simd for i in 1:problem_size
        W[i, :] = getvec(order[i], memory_effect, gen)
    end
    x1 = copy(u0)

    du = zeros(T, problem_size)
    for k in 2:N
        tk = mesh[k]
        L = Int(min(k - 1, L0))
        #if iip
        prob.f(du, x1, p, tk)
        #else
        #    du = prob.f(x1, p, tk)
        #end

        for i in 1:problem_size
            x1[i] = du[i] * halpha[i] - W[i, 2:(L + 1)]' * z[i, (k - 1):-1:(k - L)] + u0[i]
        end
        z[:, k] = x1 - u0
    end
    y = z + repeat(u0, 1, N)
    y = collect(Vector{T}, eachcol(y))

    return DiffEqBase.build_solution(prob, alg, mesh, y)
end

"""
P-th precision polynomial generate function

```math
g_p(z)=\\sum_{k=1}^p \\frac{1}{k}(1-z)^k
```
"""
function genfun(p)
    a = collect(1:(p + 1))
    A = Vandermonde(a)'
    return (1 .- a') * inv(A')
end

function getvec(order, problem_size, g)
    p = length(g) - 1
    b = 1 + order
    g0 = g[1]
    w = Float64[]
    push!(w, g[1]^order)

    @fastmath @inbounds @simd for N in 2:p
        M = N - 1
        dA = b / M
        temp = (-(g[2:N] .* collect((1 - dA):(-dA):(1 - b))))' * w[M:-1:1] / g0
        push!(w, temp)
    end

    @fastmath @inbounds @simd for k in (p + 1):problem_size
        M = k - 1
        dA = b / M
        temp = (-(g[2:(p + 1)] .* collect((1 - dA):(-dA):(1 - p * dA))))' *
               w[M:-1:(k - p)] / g0
        push!(w, temp)
    end
    return w
end
