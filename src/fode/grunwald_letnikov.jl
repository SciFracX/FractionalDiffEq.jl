@concrete struct GLCache{iip, T}
    prob
    alg
    mesh
    u0
    order
    horder
    y
    p
    kwargs
end

Base.eltype(::GLCache{iip, T}) where {iip, T} = T

function SciMLBase.__init(prob::FODEProblem, alg::GL; dt = 0.0, kwargs...)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    # GL method is only for commensurate order FODE
    (; order, u0, tspan, p) = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    order = order[1]
    horder = dt^order[1]
    n::Int64 = floor(Int64, (tfinal - t0) / dt) + 1
    iip = isinplace(prob)
    T = eltype(u0)
    mesh = collect(T, tspan[1]:dt:tspan[2])

    # Initialize solution
    result = zeros(T, length(u0), n)
    result[:, 1] .= u0

    return GLCache{iip, T}(prob, alg, mesh, u0, order, horder, result, p, kwargs)
end

function SciMLBase.solve!(cache::GLCache{iip, T}) where {iip, T}
    (; prob, alg, mesh, u0, order, horder, y, p) = cache
    prob = _is_need_convert!(prob)
    n = length(mesh)
    l = length(u0)

    # generating generalized binomial Corder
    Corder = zeros(T, n)
    Corder[1] = 1

    @fastmath @inbounds @simd for j in range(2, n, step = 1)
        Corder[j] = (1 - (1 + order) / (j - 1)) * Corder[j - 1]
    end

    du = zeros(T, l)

    @fastmath @inbounds @simd for k in range(2, n, step = 1)
        summation = zeros(T, length(u0))

        @fastmath @inbounds @simd for j in range(1, k - 1, step = 1)
            for i in eachindex(summation)
                summation[i] += Corder[j + 1] * y[i, k - j]
            end
        end

        prob.f(du, y[:, k - 1], p, mesh[k])
        y[:, k] = @. horder * du - summation
    end
    u = collect(Vector{eltype(u0)}, eachcol(y))

    return DiffEqBase.build_solution(prob, alg, mesh, u)
end
