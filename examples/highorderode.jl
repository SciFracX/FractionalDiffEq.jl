using OrdinaryDiffEq, FractionalDiffEq, Plots

function fun(du, u, p, x)
    du = u + 2 * x^2 - x - 3
end

prob = SecondOrderODEProblem(fun, 0, 0, (0, 5))
sol = OrdinaryDiffEq.solve(prob, Tsit5())

fodeprob = SingleTermFODEProblem((x, u) -> u + 2 * x^2 - x - 3, 2, [0, 0], (0, 5))
fodesol = FractionalDiffEq.solve(fodeprob, 0.001, PECE())

analytical(x) = -2 * x^2 + x + exp(-x) - 1

plot(sol, vars = (2))
plot!(fodesol, ls = :dash)

rightfun(t, x) = 0
multiprob = MultiTermsFODEProblem([1, -9, 20], [2, 1, 0], rightfun, [1, 2], (0, 1))
multisol = FractionalDiffEq.solve(multiprob, 0.01, PIPECE())
anal(x) = exp(4 * x) * (3 - 2 * exp(x))

analytics = anal.(collect(0:0.01:1))
plot(multisol)
plot!(collect(0:0.01:1), analytics, ls = :dash)
