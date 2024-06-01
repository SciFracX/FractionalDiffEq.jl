using FractionalDiffEq, Plots
α = 0.95;
ϕ(x) = 0.5;
τ = [2, 2.6]
fun(t, y, ϕ1, ϕ2) = 2 * ϕ1 / (1 + ϕ2^9.65) - y
prob = FDDEProblem(fun, ϕ, α, τ, 100)
delayed, y = solve(prob, 0.01, DelayPECE())

p1 = plot(delayed[1, :], y)
p2 = plot(delayed[2, :], y)
plot(p1, p2, layout = (1, 2))
