"""
# Usage

    solve(prob::MultiTermsFODEProblem, t, p, ClosedFormHighPercision())

Closed form high precision algorithms for multi term ordinary differential equations.
"""
struct ClosedFormHighPercision <: FractionalDiffEqAlgorithm end

function solve(prob::MultiTermsFODEProblem, h, p, ::ClosedFormHighPercision)
    @unpack parameters, orders, rightfun, rparameters, rorders, tspan =  prob
    t0 = tspan[1]; T = tspan[2]
    t = collect(t0:h:T)
    u = rightfun.(t)
    n = length(t)
    na = orders[:]
    nb = rorders[:]
    a = parameters[:]
    b = rparameters[:]
    vec = [na; nb]
    u = u[:]
    g = genfun(p)
    t = t[:]
    W = zeros(n, n)


    for i=1:length(vec)
        W[i, :] = getvec(vec[i], n, g)
    end

    D1 = b./h.^nb
    nA = length(a)
    y1 = zeros(n)

    #W=W'
    D = sum((a[:].*W[1, 1:nA])./h.^na)

    for i = 2:n
        A = y1[i-1:-1:1]'*W[2:i, 1:nA]
        y1[i] = (u[i]-sum(A.*a./h.^na))/D
    end

    y = zeros(n)

    for i = 2:n
        # Here are some problems with the indexing
        y[i] = (W[1:i, nA+1:end]*D1)'*y1[i:-1:1]
    end

    return y
end