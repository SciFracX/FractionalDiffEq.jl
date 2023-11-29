function solve(prob::FODEProblem, ::GL, dt = 0.0)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    # GL method is only for same order FODE
    @unpack f, order, u0, tspan, p = prob
    t0 = tspan[1]; T = tspan[2]
    t = collect(Float64, t0:dt:T)
    order = order[1]
    horder = dt^order[1]
    n::Int64 = floor(Int64, (T-t0)/dt)+1
    l = length(u0)

    # Initialize solution
    result = zeros(Float64, length(u0), n)
    result[:, 1] = u0

    # generating generalized binomial Corder
    Corder = zeros(Float64, n)
    Corder[1] = 1

    @fastmath @inbounds @simd for j in range(2, n, step=1)
        Corder[j] = (1-(1+order)/(j-1))*Corder[j-1]
    end

    du = zeros(Float64, l)

    @fastmath @inbounds @simd for k in range(2, n, step=1)
        summation = zeros(Float64, length(u0))

        @fastmath @inbounds @simd for j in range(1, k-1, step=1)
            for i in eachindex(summation)
                summation[i] += Corder[j+1]*result[i, k-j]
            end
        end

        f(du, result[:, k-1], p, t0+(k-1)*dt)
        result[:, k] = @. horder*du-summation
    end
    return FODESystemSolution(t, result)
end