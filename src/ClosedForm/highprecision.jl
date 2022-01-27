import FractionalDiffEq.FractionalDiffEqAlgorithm

"""
Closed form high precision algorithms for multi term ordinary differential equations
"""
struct ClosedFormHighPercision <: FractionalDiffEqAlgorithm end

function solve(prob::MultiTermsFODEProblem, t, p, ::ClosedFormHighPercision)
    a, na, rightfun, b, nb =  prob.parameters, prob.orders, prob.rightfun, prob.rparameters, prob.rorders
    u = rightfun.(t)
    h = t[2]-t[1]
    n = length(t)
    na = na[:]
    nb = nb[:]
    b = b[:]
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

function genfun(p)
    a = collect(1:p+1)
    A = Vandermonde(a)'
    return (1 .-a')*inv(A')
end

function getvec(α, n, g)
    p = length(g)-1
    b = 1 + α
    g0 = g[1]
    w = Float64[]
    push!(w, g[1]^α)

    for m = 2:p
        M = m-1
        dA = b/M
        temp = (-(g[2:m] .*collect((1-dA):-dA:(1-b))))' *w[M:-1:1]/g0
        push!(w, temp)
    end

    for k = p+1:n
        M = k-1
        dA = b/M
        temp = (-(g[2:(p+1)] .*collect((1-dA):-dA:(1-p*dA))))' *w[M:-1:(k-p)]/g0
        push!(w, temp)
    end

    return w
end