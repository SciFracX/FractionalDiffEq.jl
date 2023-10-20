using FractionalDiffEq, Plots, SpecialFunctions

phi(x)=0
α=0.9; τ = 0.1
f(t, y, ϕ)=2/gamma(3-α)*t^(2-α)-1/gamma(2-α)*t^(1-α)+2*τ*t-τ^2-τ-y+ϕ

prob=FDDEProblem(f, phi, α, τ, (2, 5))
sol=solve(prob, 0.01, DelayPI())

plot(collect(2:0.01:5), sol)

plot!(t->t^2-t, ls=:dash)