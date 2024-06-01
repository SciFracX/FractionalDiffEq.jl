using FractionalDiffEq
h = 0.01
T = 50
α = 0.97
τ = 2
q = 0.5
f(t, ϕ, y) = 2 * ϕ / (1 + ϕ^9.65) - y

prob = FDDEProblem(f, q, α, τ, T)
x, y = solve(prob, h, DelayABM())

using Plots
plot(x[500:end], y[500:end])
