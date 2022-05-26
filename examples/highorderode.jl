using OrdinaryDiffEq, FractionalDiffEq, Plots

odeprob = SecondOrderODEProblem((v, u, p, t)->1-u, 0, 0, (0, 5))
odesol = OrdinaryDiffEq.solve(odeprob, Tsit5())

fodeprob = SingleTermFODEProblem((t, u)->1-u, 2, 0, (0, 5))
fodesol = FractionalDiffEq.solve(fodeprob, 0.001, PECE())

plot(odesol)
plot!(fodesol, ls=:dash)