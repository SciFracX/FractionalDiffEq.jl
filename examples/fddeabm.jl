h=0.01
T=50
α=0.97
tau=2
ϕ=0.5

f(t, ϕ, y) = 2*ϕ/(1+ϕ^9.65)-y

x, y=solve(f, ϕ, α, τ, T, h, DelayABM())

using Plots
plot(x[500:end], y[500:end])