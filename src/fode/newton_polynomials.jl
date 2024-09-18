function solve(prob::FODEProblem, alg::NewtonPolynomial; dt = 0.0)
    dt ≤ 0 ? throw(ArgumentError("dt must be positive")) : nothing
    (; f, order, u0, tspan, p) = prob
    t0 = tspan[1]
    tfinal = tspan[2]
    order = order[1]
    t = collect(Float64, t0:dt:tfinal)
    M = 1 - order + order / gamma(order)
    N::Int = ceil(Int, (tfinal - t0) / dt)
    l = length(u0)
    result = zeros(Float64, l, N + 1)

    result[:, 1] = u0
    temp = zeros(l)
    f(temp, result[:, 1], p, t[1])
    result[:, 2] = result[:, 1] + dt .* temp
    temptemp = zeros(l)
    f(temptemp, result[:, 2], p, t[2])
    result[:, 3] = result[:, 2] + (dt / 2) .* (3 .* temptemp - temp)

    temp1 = zeros(Float64, l)
    temp2 = zeros(Float64, l)
    temp3 = zeros(Float64, l)

    for n in 3:N
        f(temp1, result[:, n], p, t[n])
        f(temp2, result[:, n - 1], p, t[n - 1])
        f(temp3, result[:, n - 2], p, t[n - 2])
        result[:, n + 1] = result[:, n] +
                           (1 - order) / M * (temp1 - temp2) +
                           order .* M .* dt .*
                           (23 / 12 * temp1 - 4 / 3 * temp2 + 5 / 12 * temp3)
    end
    u = collect(Vector{eltype(u0)}, eachcol(result))

    return DiffEqBase.build_solution(prob, alg, t, u)
end
