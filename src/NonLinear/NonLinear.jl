"""
# Usage

    solve(f, α, x0, h, t, NonLinearAlg())

Nonlinear algorithm for nonlinear fractional differential equations.

### References

Dingyu Xue, Northeastern University, China ISBN:9787030543981
"""
struct NonLinearAlg <: FractionalDiffEqAlgorithm end

function solve(prob::FODESystem, h, tn, ::NonLinearAlg, L0=1e10)
    @unpack f, α, x0 = prob
    n = length(x0)
    m = round(Int, tn/h)+1
    g = genfun(1)
    g = g[:]
    x0 = x0[:]
    ha = h.^α
    z = zeros(n, m)
    x1 = copy(x0) # Here we pass the value of x0 to x1. Honestly, I kept finding this bug for almost a whole night😅


    # All of the min(m, L0+1) is to set the memory effect.
    SetMemoryEffect = Int64(min(m, L0+1))
    W = zeros(n, SetMemoryEffect) #Initializing W a n*m matrix

    @fastmath @inbounds @simd for i = 1:n
        W[i, :] = getvec(α[i], SetMemoryEffect, g)
    end

    @fastmath @inbounds @simd for k = 2:m
        tk = (k-1)*h
        L = min(Int64(k-1), Int64(L0))
        @fastmath @inbounds @simd for i = 1:n
            x1[i] = f(tk, x1, i)*ha[i] - W[i, 2:L+1]'*z[i, k-1:-1:k-L] + x0[i]
        end
        z[:, k] = x1 - x0
    end

    result = (z + repeat(x0, 1, m))'
    
    return result

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

    @fastmath @inbounds @simd for m = 2:p
        M = m-1
        dA = b/M
        temp = (-(g[2:m] .*collect((1-dA):-dA:(1-b))))' *w[M:-1:1]/g0
        push!(w, temp)
    end

    @fastmath @inbounds @simd for k = p+1:n
        M = k-1
        dA = b/M
        temp = (-(g[2:(p+1)] .*collect((1-dA):-dA:(1-p*dA))))' *w[M:-1:(k-p)]/g0
        push!(w, temp)
    end
    return w
end