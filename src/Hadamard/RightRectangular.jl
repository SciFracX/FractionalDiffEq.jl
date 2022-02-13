struct HadamardRRect <: FractionalDiffEqAlgorithm end

function solve(f, α, h, u0, a, b, ::HadamardRRect)
    N = floor(Int, (b-a)/h)
    y = zeros(N+1)
    leftsum = zero(Float64)
    y = zeros(N+1)
    y[1]=u0

    for i=1:N-1
        for j=0:i
            leftsum += ωᵢₙ(j, N, h, α, a)*y[j+1]
        end
        # How can we handle the uₙ then?
        y[i+1]=(f(b, 0)-leftsum)/ωᵢₙ(i, N, h, α, a)
    end

    return y
end

function ωᵢₙ(i, n, h, α, x₀)
    if 0 ≤ i ≤ n-2
        return 1/gamma(1-α)*((log((x₀+n*h)/(x₀+i*h)))^(-α) - (log((x₀+n*h)/(x₀+(i+1)*h)))^(-α))
    elseif i == n-1
        return 1/gamma(1-α)*(log((x₀+n*h)/(x₀+i*h)))^(-α)
    end
end

fun(x, y) = 1-y

result=solve(fun, 0.5, 0.001, 0, 0, 1, HadamardRRect())
tspan=collect(0:0.001:1)
#analyticalfun(x)=log(x)^1.7
#analytical = analyticalfun.(tspan)

using Plots

plot(tspan, result)
#plot!(tspan, analytical)
