#https://link.springer.com/content/pdf/10.1007/s40314-019-0951-0.pdf

#struct  L1_PCM <: FractionalDiffEqAlgorithm end

function testsolve(f, α, τ, T, h, H)
    K = round(Int, τ/h)
    N = round(Int, T/h)
    #=
    un0 = 0
    for k=1:n-1
        un0 += (aₖ(n-1-k, α)-aₖ(n-k, α))*ϕ(k)
    end

    for n=1:length(-τ:h:T)

        uₙ₀ = aₖ(n-1, α)*ϕ(0) + un0

        un1 = gamma(2-α)*h^α*f(-K+n*h, un0, ϕ(-K+n*h-τ))

        un2 = gamma(2-α)*h^α*f(-K+n*h, un0+un1, ϕ(-K+n*h-τ))
    end
    =#
    result=zeros(N)
    for n=1:N
        result[n]=upredictor(result, n, α)+gamma(2-α)*h^α*f(n*h, upredictor(result, n, α)+zpredictor(result, n, f, α, K, h, H), H(-K+(n-K)*h))
    end

    return result
end

function upredictor(u, n, α)
    temp=0
    for k=1:n-1
        temp += (aₖ(n-1-k, α)-aₖ(n-k, α))*u[k]
    end
    return aₖ(n-1, α)*u[1] + temp
end

function zpredictor(u, n, f, α, K, h, H)
    return gamma(2-α)*h^α*f(-K+n*h, upredictor(u, n, α), H(-K+(n-K)*h))
end

function aₖ(k, α)
    return (k+1)^(1-α)-k^(1-α)
end
#=
f(t, u, H) = H - u + 0.2*t-0.11
function H(t)
    if t<=0
        return 0.5
    else
        return 0
    end
end

using SpecialFunctions
result=testsolve(f, 0.7, 0.1, 5, 0.01, H)
tspan = collect(0.01:0.01:5)
using Plots
plot(tspan, result)
=#