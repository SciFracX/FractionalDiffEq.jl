#https://link.springer.com/content/pdf/10.1007/s40314-019-0951-0.pdf

struct  L1_PCM <: FractionalDiffEqAlgorithm end

function testsolve(f, ϕ, α, τ, T, h, ::L1_PCM)
    K = τ/h
    N = T/h
    un0 = 0
    for k=1:n-1
        un0 += (aₖ(n-1-k, α)-aₖ(n-k, α))*ϕ(k)
    end

    for n=1:length(-τ:h:T)

        uₙ₀ = aₖ(n-1, α)*ϕ(0) + un0

        un1 = gamma(2-α)*h^α*f(-K+n*h, un0, ϕ(-K+n*h-τ))

        un2 = gamma(2-α)*h^α*f(-K+n*h, un0+un1, ϕ(-K+n*h-τ))
    end


end

function aₖ(k, α)
    return (k+1)^(1-α)-k^(1-α)
end