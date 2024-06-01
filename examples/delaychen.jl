using FractionalDiffEq, Plots
α = [0.94, 0.94, 0.94];
ϕ = [0.2, 0, 0.5];
τ = 0.009;
T = 1.4;
h = 0.001;
function delaychen!(dy, y, ϕ, t)
    a = 35
    b = 3
    c = 27
    dy[1] = a * (y[2] - ϕ[1])
    dy[2] = (c - a) * ϕ[1] - y[1] * y[3] + c * y[2]
    dy[3] = y[1] * y[2] - b * ϕ[3]
end
prob = FDDESystem(delaychen!, ϕ, α, τ, T)
#y, x=solve(prob, h, DelayABM())
#plot(x[:, 1], x[:, 2], x[:, 3], title="Fractional Order Chen Delayed System")
sol = solve(prob, h, DelayABM())
plot(sol)
