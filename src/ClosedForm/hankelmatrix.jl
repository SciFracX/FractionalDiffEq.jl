"""
# Usage

    solve(parameters, order, lparameters, lorders, u, t)

Use Closed-Form Hankel matrix algorithm to obtain numerical solution at zero initial condition.
"""
struct ClosedFormHankelM <: FractionalDiffEqAlgorithm end


function solve(prob::MultiTermsFODEProblem, t, ::ClosedFormHankelM)
    a, na, rightfun, b, nb = prob.parameters, prob.orders, prob.rightfun, prob.rparameters, prob.rorders
    h = t[2]-t[1]
    u = rightfun.(t)
    u = u[:]
    A, B = 0, 0

    g = genfun(1)
    nt = length(t)
    n = length(a)
    m = length(b)
    for i=1:n
        A = A .+ getvec(na[i], nt, g)*a[i]/(h^na[i])
    end

    for i=1:m
        B = B .+ getvec(nb[i], nt, g)*b[i]/(h^nb[i])
    end

    A = rotl90(newhankel(A[end:-1:1]))
    B = rotl90(newhankel(B[end:-1:1]))

    y = B*inv(A)*u
    return y
end


"""
P-th precision polynomial generate function

```math
g_p(z)=\\sum_{k=1}^p \\frac{1}{k}(1-z)^k
```
"""
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

function newhankel(v)
    n = length(v)
    v = v[:]

    hankelm = zeros(n, n)
    for i=1:length(v)
        hankelm[i, 1:end-i+1] = v[i:end]
    end

    return hankelm
end