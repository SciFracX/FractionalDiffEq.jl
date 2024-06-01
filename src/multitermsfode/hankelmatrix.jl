"""
# Usage

    solve(prob::MultiTermsFODEProblem, ::ClosedFormHankelM)

Use Closed-Form Hankel matrix algorithm to obtain numerical solution at zero initial condition.
"""
struct ClosedFormHankelM <: MultiTermsFODEAlgorithm end

function solve(prob::MultiTermsFODEProblem, h, ::ClosedFormHankelM)
    @unpack parameters, orders, rightfun, rparameters, rorders, u0, tspan = prob
    t0 = tspan[1]
    T = tspan[2]
    t = collect(t0:h:T)
    u = rightfun.(t)
    u = u[:]
    A, B = 0, 0

    g = genfun(1)
    nt = length(t)
    n = length(parameters)
    m = length(rparameters)
    for i in 1:n
        A = A .+ getvec(orders[i], nt, g) * parameters[i] / (h^orders[i])
    end

    for i in 1:m
        B = B .+ getvec(rorders[i], nt, g) * rparameters[i] / (h^rorders[i])
    end

    A = rotl90(newhankel(A[end:-1:1]))
    B = rotl90(newhankel(B[end:-1:1]))

    y = B * inv(A) * u
    return FODESolution(t, y)
end

function newhankel(v)
    n = length(v)
    v = v[:]

    hankelm = zeros(n, n)
    for i in 1:length(v)
        hankelm[i, 1:(end - i + 1)] = v[i:end]
    end

    return hankelm
end
